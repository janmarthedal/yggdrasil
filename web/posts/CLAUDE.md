# Posts — Finite Element Method

## Purpose

The markdown files in this folder form a series of posts introducing the finite
element method, accompanying the `yggdrasil` library in this repository.
The posts explain the mathematical theory and code constructs and reference the
accompanying library code in ./yggdrasil/ and possibly also example code in
./examples/.
*Every* function and class in the library should be referenced from at least one
post.
The posts should be prefixed with 01, 02, and so on, such that they form a
linear narrative.
Each post should be of small or medium size and as self-contained as possible.
Plenty of examples, illustrations and animations should be included.
The file `index.md` provides an introduction to and overview of the post series.

## Writing conventions

- Reference library source files by path relative to the repo root, e.g.
  `yggdrasil/assemble.py`.
- Math is written in LaTeX delimited by `$...$` (inline) and `$$...$$` (display).
- Mark incomplete sections with `**TODO**`.
- Prefer SVG over PNG for inline images when the illustration is fit for vector
  graphics.
- Consider using [three.js](https://threejs.org/) for 3D illustrations.
