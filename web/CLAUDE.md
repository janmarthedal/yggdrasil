# Posts and media on the Finite Element Method and the Yggdrasil project

## Purpose
The folder ./posts/ contains a series of posts introducing the finite element
method, accompanying the `yggdrasil` library in this repository.
The posts explain the mathematical foundation of the finite element method
and the code constructs that make up the library in relation to the theory.
There will be little abstract mathematical theory, but references to the
underlying theory should be included. Resources accessible online are preferred,
but authoritative books and papers may also be mentioned.
It is an important goal that every function and class from the library are
referenced from at least one post at some point.

## Post guidelines
- Plenty of examples, illustrations and animations should be included.
- Each post is of small size and as self-contained as possible.
- Use a shared set of notation, which is listed in `notation.md`.
  Update this file as posts are written.

## Writing conventions
- Use markdown syntax.
- Reference library source files using links, e.g., reference
  yggdrasil/assemble.py as "[`assemble.py`](LIBROOT/assemble.py)".
- Reference image files by path relative to `MEDIAROOT`.
- Reference other posts by path relative to `POSTROOT`. The post link for file
  `some-post.md` should be `POSTROOT/some-post/`.
- Math is written in LaTeX delimited by `$...$` (inline) and `$$...$$` (display).
- Prefer SVG over PNG for inline images when the illustration is fit for vector
  graphics.
- Consider using [three.js](https://threejs.org/) for 3D illustrations.
- Headings should be avoided. Let the text flow as continuous prose.

## Post Specification
The file [post specification](./post-specification.md) contains a list
of current and planned posts in chronological order.
In addition to each post's filename and title, it also includes a specification
of each post's content.
