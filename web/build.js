import markdownItKatex from '@vscode/markdown-it-katex';
import MarkdownIt from 'markdown-it';
import fs from 'node:fs';
import path from 'node:path';

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
  img {
    max-width: 100%;
    height: auto;
    display: block;
    margin: 0 auto;
  }
`;

const md = new MarkdownIt({ html: true }).use(markdownItKatex.default);

const postsDir = 'posts';
const mediaDir = 'media';
const siteDir = '_site';

function buildFile(file) {
  const src = path.join(postsDir, file);
  const basename = path.basename(file, '.md');
  const dest = basename === 'index'
    ? path.join(siteDir, 'index.html')
    : path.join(siteDir, basename, 'index.html');
  fs.mkdirSync(path.dirname(dest), { recursive: true });

  const markdown = fs.readFileSync(src, 'utf8')
    .replaceAll('MEDIAROOT', '/media')
    .replaceAll('POSTROOT', '.');
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

function copyMedia() {
  const destDir = path.join(siteDir, mediaDir);
  fs.mkdirSync(destDir, { recursive: true });
  for (const file of fs.readdirSync(mediaDir)) {
    const src = path.join(mediaDir, file);
    const dest = path.join(destDir, file);
    fs.copyFileSync(src, dest);
    console.log(`${src} → ${dest}`);
  }
}

function buildAll() {
  fs.mkdirSync(siteDir, { recursive: true });
  copyMedia();
  const files = fs.readdirSync(postsDir).filter(f => f.endsWith('.md'));
  for (const file of files) buildFile(file);
}

buildAll();

if (process.argv.includes('--watch')) {
  const { watch } = await import('chokidar');
  console.log(`Watching ${postsDir}/ and ${mediaDir}/ for changes...`);
  watch(postsDir).on('all', (event, file) => {
    const base = path.basename(file);
    if (!base.endsWith('.md')) return;
    console.log(`Post ${event}: ${file}`);
    buildFile(base);
  });
  watch(mediaDir).on('all', (event, file) => {
    console.log(`Media ${event}: ${file}`);
    copyMedia();
  });
}
