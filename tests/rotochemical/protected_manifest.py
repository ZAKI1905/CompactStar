"""Fresh manifest generation plus scratch-only topology negative controls."""

import argparse
import copy
import json
import sys
import tempfile
import subprocess
from pathlib import Path

sys.dont_write_bytecode = True
from manifest import check_manifest, generate_manifest  # noqa: E402


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    args = parser.parse_args()
    scratch_root = args.scratch_root.resolve()
    scratch_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="manifest-controls-", dir=scratch_root) as name:
        scratch = Path(name)
        entry = scratch / "entry-hashes.json"
        original = generate_manifest(args.source_root, entry)
        mutations = {
            "missing_path": lambda d: d["protected"].pop(),
            "changed_hash": lambda d: d["protected"][0].update(actual="0" * 64),
            "duplicate_path": lambda d: d["protected"].__setitem__(-1, d["protected"][0]),
            "extra_path": lambda d: d["protected"].append(copy.deepcopy(d["protected"][0])),
            "missing_baseline": lambda d: d["baselines"].pop(),
            "missing_special": lambda d: d["special"].pop(),
        }
        for label, mutate in mutations.items():
            data = copy.deepcopy(original)
            mutate(data)
            path = scratch / (label + ".json")
            path.write_text(json.dumps(data))
            try:
                check_manifest(args.source_root, path)
            except RuntimeError:
                print(label, "REFUSED")
            else:
                raise RuntimeError("manifest mutation escaped: " + label)
        check_manifest(args.source_root, entry)

        # Mutate an actual protected source only inside a disposable detached
        # worktree, after its manifest has been generated.  The current source
        # tree is never modified by this control.
        mutation_source = scratch / "mutated-source"
        subprocess.run(
            ["git", "worktree", "add", "--detach", str(mutation_source), "HEAD"],
            cwd=args.source_root, check=True, stdout=subprocess.DEVNULL,
        )
        try:
            mutation_entry = scratch / "mutation-entry.json"
            mutation = generate_manifest(mutation_source, mutation_entry)
            protected_path = mutation["protected"][0]["path"]
            with (mutation_source / protected_path).open("a") as stream:
                stream.write("\n")
            try:
                check_manifest(mutation_source, mutation_entry)
            except RuntimeError:
                print("protected_source_after_entry REFUSED")
            else:
                raise RuntimeError("protected source mutation escaped")
        finally:
            subprocess.run(
                ["git", "worktree", "remove", "--force", str(mutation_source)],
                cwd=args.source_root, check=True, stdout=subprocess.DEVNULL,
            )
    print("PASS fresh exact 33 source paths, 10 baselines, 3 special hashes")


if __name__ == "__main__":
    main()
