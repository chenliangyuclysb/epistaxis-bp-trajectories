"""Build and verify the local code-only publication archive."""

import hashlib
from importlib import import_module
from importlib.metadata import version
from pathlib import Path
import sys
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo


ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"
PACKAGE = OUTPUTS / "p6_code_release.zip"


def code_members():
    sources = [path for path in (ROOT / "code").iterdir()
               if path.suffix in (".py", ".R") and (path.stem == "config" or path.name[:2] in
                   ("01", "02", "03", "04", "05", "06", "07"))]
    paths = (*sources, ROOT / "README_code.md", ROOT / "ENVIRONMENT.md")
    expected = {}
    for path in sorted(paths):
        if path.resolve() != path or not path.is_file():
            raise ValueError("A release source is missing or not a regular in-project file.")
        relative = path.relative_to(ROOT).as_posix()
        expected[relative] = path.read_bytes()
    return expected


def write_package(members, writer):
    staged = PACKAGE.with_suffix(".zip.tmp")
    try:
        with ZipFile(staged, "w") as archive:
            for relative, content in members.items():
                item = ZipInfo(relative)
                item.compress_type = ZIP_DEFLATED
                archive.writestr(item, content, compresslevel=9)
        with ZipFile(staged) as archive:
            if archive.namelist() != list(members) or archive.testzip() is not None:
                raise ValueError("Archive structure failed the code-only check.")
            if any(archive.read(relative) != content for relative, content in members.items()):
                raise ValueError("Archive content failed verification against its source.")
        staged.replace(PACKAGE)
    finally:
        staged.unlink(missing_ok=True)
    audit = OUTPUTS / "reproducibility"
    audit.mkdir(parents=True, exist_ok=True)
    writer(audit / "p6_code_release_manifest.csv", ("archive_member", "bytes", "sha256"),
           [(relative, len(content), hashlib.sha256(content).hexdigest()) for relative, content in members.items()])


def record_packages(p0):
    p0.record_environment(None)
    library = ROOT / ".venv" / "R" / "library"
    verified = []
    for package in ("lme4", "lmerTest"):
        lines = (library / package / "DESCRIPTION").read_text(encoding="utf-8", errors="replace").splitlines()
        package_version = next(line.partition(":")[2].strip() for line in lines if line.startswith("Version:"))
        verified.append(f"{package}: {package_version}")
    with (OUTPUTS / "env.txt").open("a", encoding="utf-8") as stream:
        stream.write("\n".join((*verified, f"patsy: {version('patsy')}", f"Pillow: {version('pillow')}", "")))


def main():
    try:
        import_module("02_p1_setup").verify_project_runtime()
        p0 = import_module("01_p0_audit")
        OUTPUTS.mkdir(parents=True, exist_ok=True)
        members = code_members()
        record_packages(p0)
        write_package(members, p0.write_csv)
    except Exception as error:
        print(f"P6 code release preparation stopped ({type(error).__name__}).", file=sys.stderr)
        return 1
    print(f"Verified {len(members)} local code and documentation members.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
