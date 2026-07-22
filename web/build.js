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
  hr {
    border: none;
    text-align: center;
    overflow: visible;
  }
  hr::after {
    content: "• • •";
    display: inline-block;
    font-size: 1.2rem;
    letter-spacing: 0.5em;
    color: #999;
    background: #fff;
    padding: 0 0.5em;
  }
`;

const md = new MarkdownIt({ html: true }).use(markdownItKatex.default);

const postsDir = 'posts';
const mediaDir = 'media';
const siteDir = '_site';

function parseFrontmatter(content) {
  const match = content.match(/^---\n([\s\S]*?)\n---\n([\s\S]*)$/);
  if (!match) return { title: null, body: content };
  const meta = {};
  for (const line of match[1].split('\n')) {
    const [key, ...rest] = line.split(':');
    if (key) meta[key.trim()] = rest.join(':').trim();
  }
  return { title: meta.title ?? null, body: match[2] };
}

function buildFile(relPath) {
  const src = path.join(postsDir, relPath);
  const basename = path.basename(relPath, '.md');
  const dir = path.dirname(relPath);
  const dest = (basename === 'index' && dir === '.')
    ? path.join(siteDir, 'index.html')
    : path.join(siteDir, dir, basename, 'index.html');
  fs.mkdirSync(path.dirname(dest), { recursive: true });

  const raw = fs.readFileSync(src, 'utf8');
  const { title, body: rawBody } = parseFrontmatter(raw);
  const pageTitle = title ?? basename;

  const markdown = rawBody
    .replaceAll('MEDIAROOT', '/media')
    .replaceAll('POSTROOT', '')
    .replaceAll('LIBROOT', 'https://github.com/janmarthedal/yggdrasil/tree/main/yggdrasil');
  const bodyHtml = md.render(markdown);

  const html = `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>${pageTitle}</title>
  <link rel="stylesheet" href="${KATEX_CSS}">
  <style>${PAGE_CSS}</style>
</head>
<body>
<h1>${pageTitle}</h1>
${bodyHtml}
</body>
</html>
`;

  fs.writeFileSync(dest, html);
  console.log(`${src} → ${dest}`);
}


function getMarkdownFiles(dir, base = '') {
  const files = [];
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const rel = path.join(base, entry.name);
    if (entry.isDirectory()) {
      files.push(...getMarkdownFiles(path.join(dir, entry.name), rel));
    } else if (entry.name.endsWith('.md')) {
      files.push(rel);
    }
  }
  return files;
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
  for (const file of getMarkdownFiles(postsDir)) buildFile(file);
}

buildAll();

if (process.argv.includes('--watch')) {
  const { watch } = await import('chokidar');
  console.log(`Watching ${postsDir}/ and ${mediaDir}/ for changes...`);
  watch(postsDir).on('all', (event, file) => {
    const rel = path.relative(postsDir, file);
    if (!rel.endsWith('.md')) return;
    console.log(`Post ${event}: ${file}`);
    buildFile(rel);
  });
  watch(mediaDir).on('all', (event, file) => {
    console.log(`Media ${event}: ${file}`);
    copyMedia();
  });
}
