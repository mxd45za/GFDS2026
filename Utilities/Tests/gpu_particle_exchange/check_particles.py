"""Read actual Smokeview particle records; do not infer migration from exit status."""
import argparse
import csv
import json
import struct
from pathlib import Path


def read_particles(path):
    with path.open('rb') as stream:
        def record():
            head = stream.read(4)
            if not head:
                return None
            size, = struct.unpack('<i', head)
            data = stream.read(size)
            tail = stream.read(4)
            assert len(data) == size and tail == head, f'Truncated record: {path}'
            return data

        def integers(data):
            return struct.unpack(f'<{len(data)//4}i', data)

        def floats(data):
            return struct.unpack(f'<{len(data)//4}f', data)

        assert integers(record()) == (1,), 'Unsupported endian format'
        record()  # Version
        nclasses, = integers(record())
        labels = []
        for _ in range(nclasses):
            nq, _ = integers(record())
            names = []
            for _ in range(nq):
                names.append(record().decode().strip())
                record()  # Units
            labels.append(names)
        frames = []
        while (data := record()) is not None:
            time, = floats(data)
            particles = {}
            for class_id, names in enumerate(labels, 1):
                count, = integers(record())
                xyz = floats(record())
                tags = integers(record())
                values = floats(record()) if names else ()
                assert len(xyz) == 3*count and len(tags) == count
                for i, tag in enumerate(tags):
                    assert tag not in particles, f'Duplicate tag {tag} in {path}'
                    particles[tag] = dict(class_id=class_id,
                        xyz=[xyz[i + j*count] for j in range(3)],
                        values={name: values[i+j*count] for j, name in enumerate(names)})
            frames.append((time, particles))
        return frames


def check_case(folder, chid):
    log = (folder / 'run.log').read_text()
    assert 'FDS completed successfully' in log and 'ERROR' not in log, f'Failed run: {folder}'
    paths = sorted(folder.glob(f'{chid}_*.prt5'))
    assert paths, f'No particle output in {folder}'
    meshes = [read_particles(path) for path in paths]
    assert len({len(mesh) for mesh in meshes}) == 1
    assert len(meshes[0]) > 2
    known, owners, crossed = {}, {}, set()
    max_conservation_error = 0.0
    history = []
    for frame in zip(*meshes):
        time = frame[0][0]
        assert all(abs(time-f[0]) < 1e-6 for f in frame)
        current = {}
        for owner, (_, particles) in enumerate(frame, 1):
            for tag, p in particles.items():
                assert tag not in current, f'Duplicate tag {tag} at t={time}'
                xlow = (owner-1) % 2 if chid == 'multidest' else owner-1
                ylow = (owner-1) // 2 if chid == 'multidest' else 0
                assert xlow-1e-5 <= p['xyz'][0] <= xlow+1+1e-5, (tag, owner, p)
                assert ylow-1e-5 <= p['xyz'][1] <= ylow+1+1e-5, (tag, owner, p)
                assert all(abs(v) < float('inf') for v in p['xyz']+list(p['values'].values()))
                if tag not in known:
                    known[tag], owners[tag] = p, owner
                elif owners[tag] != owner:
                    crossed.add(tag)
                assert p['class_id'] == known[tag]['class_id'], f'Class changed for {tag}'
                for name, value in p['values'].items():
                    if 'mass' in name.lower() or 'weight' in name.lower():
                        base = known[tag]['values'][name]
                        err = abs(value-base)/max(abs(base), 1e-30)
                        max_conservation_error = max(max_conservation_error, err)
                        assert err < 2e-6, f'{name} changed for {tag}: {err}'
                current[tag] = p
        if chid != 'removal':
            assert current.keys() == known.keys(), f'Tag lost at t={time}'
        history.append(dict(time=time, counts=[len(f[1]) for f in frame]))
    final_counts = history[-1]['counts']
    if chid in ('one_way', 'linked'):
        count = {'one_way': 32, 'linked': 1}[chid]
        assert len(known) == count and len(crossed) == count and final_counts == [0, count]
        if chid == 'linked':
            with (folder / f'{chid}_devc.csv').open() as stream:
                rows = list(csv.reader(stream))
            assert float(rows[-1][1]) > 1.3, 'INIT-linked DEVC lost the migrated particle'
    elif chid == 'many_particles':
        assert len(known) == 649 and len(crossed) == 512 and final_counts == [137, 512]
    elif chid == 'bidirectional':
        assert len(known) == 86 and len(crossed) == 86 and final_counts == [21, 43, 22]
    elif chid == 'multidest':
        assert len(known) == 123 and len(crossed) == 123 and final_counts == [0, 40, 41, 42]
    elif chid == 'periodic':
        assert len(known) == 32 and len(crossed) == 32 and final_counts == [32, 0]
    elif chid == 'removal':
        assert len(known) == 16 and sum(final_counts) == 0
        assert len(crossed) == 8
        assert all(sum(f['counts']) == 0 for f in history if f['time'] > 0.045)
    elif chid in ('reinsertion', 'count_devices'):
        assert len(known) >= 24 and final_counts[1] >= 24
        assert sum(final_counts) == len(known)
        assert all(sum(history[i]['counts']) <= sum(history[i+1]['counts']) for i in range(len(history)-1))
    else:
        raise ValueError(chid)
    if chid == 'count_devices':
        with (folder / f'{chid}_devc.csv').open() as stream:
            rows = list(csv.reader(stream))
        samples = [[float(value) for value in row] for row in rows[2:]]
        assert len(samples) > len(history), 'Counts must be sampled between particle output frames'
        assert all(value >= 0 and value.is_integer() for row in samples for value in row[1:])
        assert samples[-1][1:] == final_counts
        for frame in history:
            row = min(samples, key=lambda row: abs(row[0]-frame['time']))
            assert abs(row[0]-frame['time']) < 1e-6
            assert row[1:] == frame['counts'], (row, frame)
    return dict(case=chid, folder=str(folder), tags=len(known), observed_crossings=len(crossed),
                final_counts=final_counts, max_mass_or_weight_relative_error=max_conservation_error,
                history=history)


def compare_reference(folder, reference, chid):
    maximum_position = 0.0
    maximum_fields = {}
    for path in sorted(reference.glob(f'{chid}_*.prt5')):
        baseline = read_particles(path)
        result = read_particles(folder / path.name)
        assert len(baseline) == len(result)
        for (t, a), (u, b) in zip(baseline, result):
            assert abs(t-u) < 1e-6 and a.keys() == b.keys(), (chid, path.name, t, len(a), len(b))
            for tag in a:
                maximum_position = max(maximum_position,
                    *(abs(x-y) for x, y in zip(a[tag]['xyz'], b[tag]['xyz'])))
                for name, value in a[tag]['values'].items():
                    error = abs(value-b[tag]['values'][name])
                    maximum_fields[name] = max(maximum_fields.get(name, 0.0), error)
    assert maximum_position < 5e-5, maximum_position
    assert all(error < 1e-4 for error in maximum_fields.values()), maximum_fields
    if chid == 'count_devices':
        def device_rows(root):
            with (root / f'{chid}_devc.csv').open() as stream:
                rows = list(csv.reader(stream))
            return [[float(value) for value in row] for row in rows[2:]]
        baseline = device_rows(reference)
        result = device_rows(folder)
        assert len(baseline) == len(result)
        for a, b in zip(baseline, result):
            assert abs(a[0]-b[0]) < 1e-6 and a[1:] == b[1:], (a, b)
    return dict(reference=str(reference), maximum_position_m=maximum_position,
                maximum_field_absolute_errors=maximum_fields)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('folder', type=Path)
    parser.add_argument('--case', default='one_way')
    parser.add_argument('--reference', type=Path)
    args = parser.parse_args()
    result = check_case(args.folder, args.case)
    if args.reference:
        result['comparison'] = compare_reference(args.folder, args.reference, args.case)
    print(json.dumps(result, indent=2))
