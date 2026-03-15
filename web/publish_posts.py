#!/usr/bin/env python3
"""Copy web/posts/ and web/media/ to target directories, rewriting placeholder paths."""

import argparse
import shutil
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
POSTS_SRC = SCRIPT_DIR / "posts"
MEDIA_SRC = SCRIPT_DIR / "media"


def confirm_overwrite(path: Path) -> bool:
    response = input(f"Overwrite {path}? [y/N] ").strip().lower()
    return response == "y"


def copy_file(src: Path, dst: Path, force: bool) -> bool:
    if dst.exists() and not force:
        if not confirm_overwrite(dst):
            return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def rewrite_placeholders(path: Path, replacements: dict[str, str]) -> None:
    text = path.read_text(encoding="utf-8")
    for placeholder, value in replacements.items():
        text = text.replace(placeholder, value)
    path.write_text(text, encoding="utf-8")


def publish(
    posts_dest: Path,
    media_dest: Path,
    postroot: str,
    mediaroot: str,
    libroot: str | None,
    force: bool,
) -> None:
    replacements = {
        "POSTROOT": postroot,
        "MEDIAROOT": mediaroot,
    }
    if libroot is not None:
        replacements["LIBROOT"] = libroot

    # Copy media files as-is
    for src in MEDIA_SRC.rglob("*"):
        if src.is_file():
            rel = src.relative_to(MEDIA_SRC)
            dst = media_dest / rel
            copy_file(src, dst, force)

    # Copy post files and rewrite placeholders in markdown files
    for src in POSTS_SRC.rglob("*"):
        if src.is_file():
            rel = src.relative_to(POSTS_SRC)
            dst = posts_dest / rel
            if copy_file(src, dst, force) and src.suffix == ".md":
                rewrite_placeholders(dst, replacements)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Publish posts and media to a target repository."
    )
    parser.add_argument(
        "--posts-dest",
        required=True,
        type=Path,
        help="Destination directory for posts",
    )
    parser.add_argument(
        "--media-dest",
        required=True,
        type=Path,
        help="Destination directory for media files",
    )
    parser.add_argument(
        "--postroot",
        required=True,
        help="Replacement value for POSTROOT placeholder",
    )
    parser.add_argument(
        "--mediaroot",
        required=True,
        help="Replacement value for MEDIAROOT placeholder",
    )
    parser.add_argument(
        "--libroot",
        default=None,
        help="Replacement value for LIBROOT placeholder (e.g. a GitHub URL)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing files without confirmation",
    )
    args = parser.parse_args()

    if not POSTS_SRC.is_dir():
        print(f"Error: posts source directory not found: {POSTS_SRC}", file=sys.stderr)
        sys.exit(1)
    if not MEDIA_SRC.is_dir():
        print(f"Error: media source directory not found: {MEDIA_SRC}", file=sys.stderr)
        sys.exit(1)

    publish(
        posts_dest=args.posts_dest,
        media_dest=args.media_dest,
        postroot=args.postroot,
        mediaroot=args.mediaroot,
        libroot=args.libroot,
        force=args.force,
    )


if __name__ == "__main__":
    main()
