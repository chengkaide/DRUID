/* ══════════════════════════════════════════════════════════════════════════
   druid Web 界面交互逻辑（原生 JS，无依赖）
   --------------------------------------------------------------------------
   三个状态变量撑起整个页面：
     curDir   当前选中的批次目录
     curTid   当前任务 ID（轮询用）
     since    日志游标 —— 每次轮询只取新行，不全量重传
   ══════════════════════════════════════════════════════════════════════════ */

let curDir = '';
let curTid = null;
let since = 0;
let pollTimer = null;
let tStart = 0;

// ── 小工具 ────────────────────────────────────────────────────────────────
const $ = (id) => document.getElementById(id);

/** 转义 HTML：日志里可能含 < > 等字符，不转义会把页面结构打乱 */
const esc = (s) => String(s).replace(/[&<>"']/g, (c) =>
  ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

async function api(path, opts) {
  const r = await fetch(path, opts);
  if (!r.ok) throw new Error(`HTTP ${r.status}`);
  return r.json();
}

const post = (path, body) => api(path, {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify(body)
});

// ── 页面就绪 ──────────────────────────────────────────────────────────────
document.addEventListener('DOMContentLoaded', () => {
  $('srvAddr').textContent = `本地服务 ${location.host}（仅本机可访问）`;

  $('btnBrowseDir').onclick = () => openDirBrowser(curDir || '');
  $('btnBrowseDirect').onclick = pastePath;
  $('btnCloseDir').onclick = () => $('mask').classList.remove('show');
  $('mask').onclick = (e) => { if (e.target === $('mask')) $('mask').classList.remove('show'); };
  $('btnUp').onclick = () => { if (dirCache.parent) loadDir(dirCache.parent); };
  $('btnDrives').onclick = () => loadDir('');
  $('btnPickThis').onclick = () => {
    $('mask').classList.remove('show');
    if (curPathOfBrowser) chooseDir(curPathOfBrowser);
  };
  $('btnInspect').onclick = () => chooseDir(curDir);
  $('btnRun').onclick = startRun;
  setupDrop();
});

// ══════════════════════════════════════════════════════════════════════════
// 一、选目录
// ══════════════════════════════════════════════════════════════════════════
let dirCache = {};
let curPathOfBrowser = '';

function pastePath() {
  const p = prompt('输入批次目录的完整路径，例如：\nD:\\data\\EX2022A', curDir);
  if (p && p.trim()) chooseDir(p.trim());
}

function openDirBrowser(start) {
  $('mask').classList.add('show');
  loadDir(start);
}

async function loadDir(path) {
  $('dlist').innerHTML = '<div class="drow"><span class="nm">读取中…</span></div>';
  let info;
  try {
    info = await api('/api/browse?path=' + encodeURIComponent(path));
  } catch (e) {
    $('dlist').innerHTML = `<div class="drow"><span class="nm">读取失败：${esc(e.message)}</span></div>`;
    return;
  }
  dirCache = info;
  curPathOfBrowser = info.cwd || '';
  $('crumb').textContent = info.cwd ? `当前：${info.cwd}` : '磁盘根目录';

  const rows = [];
  if (info.parent) {
    rows.push(`<div class="drow" data-go="${esc(info.parent)}">
        <span class="nm">↩ 上一级</span><span class="tg"></span></div>`);
  }
  for (const d of (info.dirs || [])) {
    const isBatch = d.has_list && d.n_csv > 0;
    rows.push(`<div class="drow ${isBatch ? 'batch' : ''}" data-go="${esc(d.path)}">
        <span class="nm">${esc(d.name)}</span>
        <span class="tg">${isBatch ? `批次 · ${d.n_csv} CSV` : (d.n_csv ? `${d.n_csv} CSV` : '')}</span></div>`);
  }
  if (!rows.length) rows.push('<div class="drow"><span class="nm">（没有子文件夹）</span></div>');
  $('dlist').innerHTML = rows.join('');

  $('dlist').querySelectorAll('[data-go]').forEach((el) => {
    el.onclick = () => loadDir(el.getAttribute('data-go'));
  });
}

async function chooseDir(path) {
  curDir = path;
  $('chosenBox').hidden = false;
  $('btnInspect').hidden = false;
  $('chosenPath').textContent = path;
  $('facts').innerHTML = '<span class="spin"></span> 正在检查…';
  $('problems').hidden = true;
  setStepEnabled(2, false);
  setStepEnabled(3, false);

  let info;
  try {
    info = await post('/api/inspect', { dir: path });
  } catch (e) {
    say([`检查失败：${e.message}`], true);
    return;
  }

  // ── 事实标签 ──
  const f = [];
  f.push(`<span class="fact ${info.list_exists ? '' : 'bad'}">${info.list_exists ? '✓' : '✗'} 序列表</span>`);
  if (info.csv_count) f.push(`<span class="fact">✓ <b>${info.csv_count}</b> 个数据文件</span>`);
  else f.push(`<span class="fact bad">✗ 没有 CSV 数据文件</span>`);

  if (info.sequence) {
    for (const [k, v] of Object.entries(info.sequence.counts)) {
      f.push(`<span class="fact">${esc(k)} <b>${v}</b></span>`);
    }
  }
  $('facts').innerHTML = f.join('');

  // ── 问题清单 ──
  if (info.problems && info.problems.length) {
    $('problems').hidden = false;
    $('problemList').innerHTML = info.problems.map((p) => `<li>${esc(p)}</li>`).join('');
  }

  const ok = info.list_exists && info.csv_count > 0;
  setStepEnabled(2, ok);
  setStepEnabled(3, ok);
  $('btnRun').disabled = !ok;
  $('runNote').textContent = ok
    ? `将输出到：${info.out_excel}`
    : '上面列出的问题需要先解决，否则无法处理。';

  // 非阻断提示（可能有多条：缺主标/缺监控标样/剥蚀不均/仅标样…）
  const notices = info.notices && info.notices.length ? info.notices
    : (info.notice ? [info.notice] : []);
  if (ok && notices.length) {
    $('runNote').innerHTML = notices.map((n) => `• ${esc(n)}`).join('<br>');
    say(notices);
  }

  if (ok) say([`已选择：${path}`, `共 ${info.sequence ? info.sequence.rows : info.csv_count} 个测点待处理。`]);
}

function setStepEnabled(n, on) {
  const el = $('s' + n);
  if (el) el.classList.toggle('dim', !on);
}

// ══════════════════════════════════════════════════════════════════════════
// 二、拖拽上传
// ══════════════════════════════════════════════════════════════════════════
function setupDrop() {
  const z = $('dropZone');
  const stop = (e) => { e.preventDefault(); e.stopPropagation(); };
  ['dragenter', 'dragover'].forEach((ev) => z.addEventListener(ev, (e) => { stop(e); z.classList.add('hot'); }));
  ['dragleave', 'drop'].forEach((ev) => z.addEventListener(ev, (e) => { stop(e); z.classList.remove('hot'); }));

  z.onclick = async () => {
    // 弹一个隐藏的 file input，让习惯点按的用户也能用
    const inp = document.createElement('input');
    inp.type = 'file'; inp.multiple = true;
    inp.accept = '.csv,.xls,.xlsx';
    inp.onchange = () => doUpload([...inp.files]);
    inp.click();
  };

  z.addEventListener('drop', (e) => {
    const fs = [...(e.dataTransfer.files || [])];
    if (fs.length) doUpload(fs);
  });
}

async function doUpload(files) {
  $('dropZone').textContent = `正在上传 ${files.length} 个文件…`;
  const packed = [];
  for (const f of files) {
    const buf = await f.arrayBuffer();
    // base64：避开手写 multipart 解析，也避开各浏览器编码差异
    let bin = '';
    const bytes = new Uint8Array(buf);
    for (let i = 0; i < bytes.length; i++) bin += String.fromCharCode(bytes[i]);
    packed.push({ name: f.name, b64: btoa(bin) });
  }

  try {
    // 用第一个文件名推断批次名：EX2022A_12.csv → EX2022A
    const guess = (files[0].name.match(/^(.+?)_\d+\.csv$/i) || [null, '上传批次'])[1];
    const r = await post('/api/upload', { files: packed, batch_name: guess });
    $('dropZone').textContent = `也可以把 CSV / LIST 文件拖到这里上传（会临时拼成一个批次目录）`;
    if (r.ok) chooseDir(r.dir);
    else say([`上传失败：${r.error || '未知错误'}`], true);
  } catch (e) {
    say([`上传失败：${e.message}`], true);
  }
}

// ══════════════════════════════════════════════════════════════════════════
// 三、跑任务
// ══════════════════════════════════════════════════════════════════════════
async function startRun() {
  const num = (id, d) => { const v = parseFloat($(id).value); return isNaN(v) ? d : v; };
  const payload = {
    data_dir: curDir,
    bulk: $('pBulk').value,
    primary: $('pPrimary').value.trim() || '91500',
    secondary: $('pSecondary').value.trim() || 'Ple',
    plot: $('pPlot').checked,
    do_depth: $('pDepth').checked,
    n_sigma_common_pb: num('pNSigma', 2.0),
    win: num('pWin', 4.0),
    step: num('pStep', 1.0),
    trim: num('pTrim', 1.5),
    blank_dur: num('pBlank', 15.0),
    deadtime_ns: num('pDead', 0),
  };
  const out = $('pOut').value.trim();
  if (out) payload.out_excel = out;

  $('btnRun').disabled = true;
  $('term').innerHTML = '';
  since = 0;
  $('outputs').innerHTML = '';
  $('stats').hidden = true;
  $('outNote').textContent = '正在处理…';
  tStart = Date.now();

  let r;
  try {
    r = await post('/api/run', payload);
  } catch (e) {
    say([`启动失败：${e.message}`], true);
    $('btnRun').disabled = false;
    return;
  }
  if (!r.ok) { say([`启动失败：${r.error}`], true); $('btnRun').disabled = false; return; }

  curTid = r.tid;
  $('logState').innerHTML = '<span class="spin"></span> 运行中';
  pollTimer = setInterval(poll, 700);   // 700ms：肉眼足够连贯，又不至于刷爆
  poll();
}

async function poll() {
  if (!curTid) return;
  let s;
  try {
    s = await api(`/api/status?tid=${curTid}&since=${since}`);
  } catch (e) {
    return;   // 网络抖动就等下一拍，不要因此中断轮询
  }

  // 追加日志
  if (s.new_lines && s.new_lines.length) {
    say(s.new_lines);
    since += s.new_lines.length;
  }

  $('barFill').style.width = s.progress + '%';
  $('pctTxt').textContent = s.progress + '%';
  const sec = ((Date.now() - tStart) / 1000);
  $('clockTxt').textContent = sec.toFixed(0) + ' s';

  if (s.state === 'done' || s.state === 'error') {
    clearInterval(pollTimer);
    pollTimer = null;
    $('barFill').classList.add(s.state === 'done' ? 'done' : 'err');
    $('logState').textContent = s.state === 'done' ? '已完成' : '失败';
    $('btnRun').disabled = false;
    renderOutputs(s);
    if (s.state === 'error') say([s.error || '未知错误'], true);
  }
}

// ══════════════════════════════════════════════════════════════════════════
// 四、日志渲染
// ══════════════════════════════════════════════════════════════════════════
function say(lines, isBad) {
  const t = $('term');
  if (t.querySelector('.empty')) t.innerHTML = '';   // 清掉占位提示
  const frag = lines.map((ln) => {
    let cls = 'line';
    if (isBad) cls += ' bad';
    else if (/^\s*\[\d\]/.test(ln)) cls += ' stage';
    else if (/完成|结果已写出|✓/.test(ln)) cls += ' good';
    else if (/!!|错误|失败|Error|Traceback/.test(ln)) cls += ' bad';
    return `<span class="${cls}">${esc(ln)}</span>`;
  }).join('\n');
  t.appendChild(document.createRange().createContextualFragment(frag));
  t.scrollTop = t.scrollHeight;   // 始终贴到底部，像 tail -f
}

// ══════════════════════════════════════════════════════════════════════════
// 五、产物与摘要
// ══════════════════════════════════════════════════════════════════════════
function renderOutputs(s) {
  setStepEnabled(5, true);

  // ── 摘要卡片 ──
  const sm = s.summary || {};
  if (Object.keys(sm).length) {
    const cards = [];
    if (sm.n_samples) cards.push(['样品测点', sm.n_samples, '']);
    if (sm.age_median !== undefined) {
      cards.push(['年龄中位 (Ma)', sm.age_median, '']);
      cards.push(['5–95% 区间', `${sm.age_p5} ~ ${sm.age_p95}`, '']);
    }
    if (sm.concordance_median !== undefined) {
      cards.push(['协和度中位', sm.concordance_median + '%',
                  sm.concordance_in_range >= 85 ? 'ok' : 'bad']);
    }
    if (sm.multi_domain !== undefined && sm.n_samples) {
      cards.push(['多年龄域样品', `${sm.multi_domain} / ${sm.n_samples}`, '']);
    }
    $('stats').hidden = false;
    $('stats').innerHTML = cards.map(([k, v, c]) =>
      `<div class="stat"><div class="k">${esc(k)}</div><div class="v ${c}">${esc(v)}</div></div>`).join('');

    // ── QC 表：标样偏差是判断这批数据能不能用的第一道关 ──
    let qcHtml = '';
    for (const q of (sm.qc || [])) {
      const dev = q['偏差_pct'] !== undefined && q['偏差_pct'] !== null ? parseFloat(q['偏差_pct']) : null;
      const cls = dev === null ? '' : (Math.abs(dev) < 3 ? 'ok' : 'bad');
      qcHtml += `<div class="stat"><div class="k">${esc(q['标样'] || '标样')} QC</div>
        <div class="v ${cls}">${dev === null ? '—' : dev.toFixed(2) + '%'}</div></div>`;
    }
    if (qcHtml) $('stats').insertAdjacentHTML('beforeend', qcHtml);
  }

  // ── 产物清单 ──
  const box = $('outputs');
  box.innerHTML = '';
  (s.outputs || []).forEach((o, i) => {
    const row = document.createElement('div');
    row.className = 'out';
    row.innerHTML = `
      <span class="ico">${o.kind === 'pdf' ? '▲' : '▤'}</span>
      <div class="meta">
        <div class="n">${esc(o.label)}</div>
        <div class="p">${esc(o.path)}</div>
      </div>
      <div class="gacts">
        <button data-dl="${i}">下载</button>
        <button data-rv="${esc(o.path)}">在文件夹中显示</button>
      </div>`;
    box.appendChild(row);
  });
  box.querySelectorAll('[data-dl]').forEach((b) => {
    b.onclick = () => { location.href = `/api/download?tid=${curTid}&i=${b.getAttribute('data-dl')}`; };
  });
  box.querySelectorAll('[data-rv]').forEach((b) => {
    b.onclick = async () => {
      const r = await post('/api/reveal', { path: b.getAttribute('data-rv') }).catch(() => null);
      if (r && !r.ok) alert('打不开：' + (r.error || '未知错误'));
    };
  });

  $('outNote').textContent = (s.outputs || []).length
    ? '点「在文件夹中显示」会直接打开资源管理器并选中文件。'
    : '本次运行没有产生可下载的文件。';
}
