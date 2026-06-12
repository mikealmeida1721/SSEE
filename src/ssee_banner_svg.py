"""
Animated SVG banner for the README (GitHub-safe: CSS animations embedded in SVG,
served as <img> — no JavaScript needed).

Theme: "perfection -> imperfection". A clean geometric phi/pi lattice on the left
gradually deforms into organic, jittering glyphs on the right; a golden
logarithmic spiral (growth factor phi per quarter turn) draws itself over a
twinkling starfield.

Output: assets/banner_ssee.svg
"""
import math
import random
import os

random.seed(42)

W, H = 1200, 300
PHI = (1 + 5**0.5) / 2

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'assets')
os.makedirs(OUT_DIR, exist_ok=True)
OUT = os.path.join(OUT_DIR, 'banner_ssee.svg')

# ── Starfield ────────────────────────────────────────────────────────────────
stars = []
for i in range(90):
    x = random.uniform(0, W)
    y = random.uniform(0, H)
    r = random.choice([0.7, 0.9, 1.1, 1.4])
    cls = f"tw{random.randint(1, 3)}"
    delay = random.uniform(0, 4)
    stars.append(
        f'<circle cx="{x:.0f}" cy="{y:.0f}" r="{r}" class="{cls}" '
        f'style="animation-delay:{delay:.1f}s"/>'
    )

# ── Golden logarithmic spiral (r = a·φ^(2θ/π)), self-drawing ────────────────
cx, cy = 985, 158
a = 2.2
pts = []
theta = 0.0
while True:
    r = a * PHI ** (2 * theta / math.pi)
    if r > 130:
        break
    pts.append((cx + r * math.cos(theta), cy + r * math.sin(theta)))
    theta += 0.05
spiral_d = 'M ' + ' L '.join(f'{x:.1f},{y:.1f}' for x, y in pts)
# total length for the dash animation
length = sum(
    math.hypot(pts[i + 1][0] - pts[i][0], pts[i + 1][1] - pts[i][1])
    for i in range(len(pts) - 1)
)
L = int(length) + 2

# ── phi/pi glyph lattice: perfect (left) -> imperfect (right) ───────────────
glyphs = []
keyframes = []
cols, rows = 16, 5
x0, dx = 60, 68
y0, dy = 58, 50
for c in range(cols):
    t = c / (cols - 1)                  # 0 = perfect, 1 = imperfect
    amp = (t ** 2) * 7.0                # jitter amplitude in px
    rot = (t ** 2) * 9.0                # max rotation in deg
    dur = 5.5 - 2.5 * t                 # faster chaos to the right
    # color: pure gold -> organic teal
    g = (int(245 - 120 * t), int(197 - 10 * t), int(66 + 120 * t))
    color = f'rgb({g[0]},{g[1]},{g[2]})'
    op = 0.16 + 0.10 * t
    if amp > 0.15:
        keyframes.append(
            f'@keyframes drift{c}{{'
            f'0%,100%{{transform:translate(0px,0px) rotate(0deg)}}'
            f'25%{{transform:translate({amp:.1f}px,{-amp*0.7:.1f}px) rotate({rot:.1f}deg)}}'
            f'55%{{transform:translate({-amp*0.8:.1f}px,{amp*0.5:.1f}px) rotate({-rot*0.7:.1f}deg)}}'
            f'80%{{transform:translate({amp*0.4:.1f}px,{amp*0.9:.1f}px) rotate({rot*0.5:.1f}deg)}}'
            f'}}'
        )
        anim = f'animation:drift{c} {dur:.1f}s ease-in-out infinite;'
    else:
        anim = ''
    for rrow in range(rows):
        ch = 'φ' if (c + rrow) % 2 == 0 else 'π'
        x = x0 + c * dx
        y = y0 + rrow * dy
        delay = (c * 0.13 + rrow * 0.31) % 2.0
        glyphs.append(
            f'<text x="{x}" y="{y}" class="gl" '
            f'style="fill:{color};opacity:{op:.2f};{anim}'
            f'animation-delay:{delay:.2f}s" '
            f'transform-origin="{x}px {y}px">{ch}</text>'
        )

svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="Georgia, 'Times New Roman', serif">
<defs>
  <radialGradient id="bg" cx="50%" cy="40%" r="90%">
    <stop offset="0%" stop-color="#0d1033"/>
    <stop offset="55%" stop-color="#070a22"/>
    <stop offset="100%" stop-color="#03040f"/>
  </radialGradient>
  <linearGradient id="gold" x1="0" y1="0" x2="1" y2="0">
    <stop offset="0%" stop-color="#f5c542"/>
    <stop offset="100%" stop-color="#ffe9a8"/>
  </linearGradient>
  <filter id="glow" x="-40%" y="-40%" width="180%" height="180%">
    <feGaussianBlur stdDeviation="3.5" result="b"/>
    <feMerge><feMergeNode in="b"/><feMergeNode in="SourceGraphic"/></feMerge>
  </filter>
</defs>
<style>
  .gl {{ font-size: 30px; }}
  circle {{ fill: #ffffff; }}
  .tw1 {{ animation: tw 2.6s ease-in-out infinite; }}
  .tw2 {{ animation: tw 3.4s ease-in-out infinite; }}
  .tw3 {{ animation: tw 4.5s ease-in-out infinite; }}
  @keyframes tw {{ 0%,100% {{opacity:.12}} 50% {{opacity:.85}} }}
  .spiral {{
    fill: none; stroke: url(#gold); stroke-width: 2.2; stroke-linecap: round;
    stroke-dasharray: {L}; opacity:.85;
    animation: draw 11s ease-in-out infinite;
  }}
  @keyframes draw {{
    0%   {{ stroke-dashoffset: {L}; opacity:.0; }}
    8%   {{ opacity:.85; }}
    62%  {{ stroke-dashoffset: 0; opacity:.85; }}
    88%  {{ stroke-dashoffset: 0; opacity:.85; }}
    100% {{ stroke-dashoffset: 0; opacity:0; }}
  }}
  {''.join(keyframes)}
  .title  {{ fill:url(#gold); font-size:54px; font-weight:bold; letter-spacing:2px; }}
  .sub    {{ fill:#cdd6f4; font-size:19px; letter-spacing:1px; }}
  .nums   {{ fill:#8fa3c8; font-size:14.5px; font-family:'Courier New',monospace; }}
  .motto  {{ fill:#6b7aa0; font-size:13px; font-style:italic; }}
</style>

<rect width="{W}" height="{H}" fill="url(#bg)"/>
{''.join(stars)}
{''.join(glyphs)}
<path class="spiral" d="{spiral_d}"/>

<g text-anchor="middle">
  <text x="540" y="128" class="title" filter="url(#glow)">SSEE — V3.6</text>
  <text x="540" y="166" class="sub">Structural Self-Energy Expansion · dark energy from φ and π</text>
  <text x="540" y="200" class="nums">w₀ = −0.840 · wₐ = −0.670 · ~3 effective parameters · 0.05σ vs DESI DR2</text>
  <text x="540" y="262" class="motto">geometric perfection underlies it — imperfection lets matter evolve</text>
</g>
</svg>
'''

with open(OUT, 'w') as f:
    f.write(svg)
print(f"Saved: {OUT}  ({os.path.getsize(OUT)/1024:.1f} KB, spiral length {L}px)")
