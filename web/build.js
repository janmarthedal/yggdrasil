import fs from 'fs';
import path from 'path';
import MarkdownIt from 'markdown-it';
import markdownItKatex from '@vscode/markdown-it-katex';

const KATEX_CSS = 'https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css';

const PAGE_CSS = `
  body {
    font-family: Georgia, serif;
    max-width: 800px;
    margin: 2rem auto;
    padding: 0 1rem;
    line-height: 1.6;
    color: #222;
  }
  code {
    font-family: monospace;
    background: #f4f4f4;
    padding: 0.1em 0.3em;
    border-radius: 3px;
  }
  pre code {
    display: block;
    padding: 1em;
    overflow-x: auto;
  }
`;

const md = new MarkdownIt({ html: true }).use(markdownItKatex.default);

const postsDir = 'posts';
const siteDir = '_site';

function buildFile(file) {
  const src = path.join(postsDir, file);
  const basename = path.basename(file, '.md');
  const dest = basename === 'index'
    ? path.join(siteDir, 'index.html')
    : path.join(siteDir, basename, 'index.html');
  fs.mkdirSync(path.dirname(dest), { recursive: true });

  const markdown = fs.readFileSync(src, 'utf8');
  const body = md.render(markdown);

  const html = `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>${basename}</title>
  <link rel="stylesheet" href="${KATEX_CSS}">
  <style>${PAGE_CSS}</style>
</head>
<body>
${body}
</body>
</html>
`;

  fs.writeFileSync(dest, html);
  console.log(`${src} → ${dest}`);
}

function buildAll() {
  fs.mkdirSync(siteDir, { recursive: true });
  const files = fs.readdirSync(postsDir).filter(f => f.endsWith('.md') && f !== 'CLAUDE.md');
  for (const file of files) buildFile(file);
}

buildAll();

if (process.argv.includes('--watch')) {
  const { watch } = await import('chokidar');
  console.log(`Watching ${postsDir}/ for changes...`);
  watch(postsDir).on('change', file => {
    const base = path.basename(file);
    if (!base.endsWith('.md') || base === 'CLAUDE.md') return;
    console.log(`Changed: ${file}`);
    buildFile(base);
  });
}
