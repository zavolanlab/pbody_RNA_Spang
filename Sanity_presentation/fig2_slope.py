import subprocess, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from PIL import Image
PDF = '/scicore/home/zavolan/GROUP/pbody_RNA_spang_zavolan_v2/materials/s41587-021-00875-x.pdf'
DPI = 600
X0, X1, Y0, Y1 = 80, 525, 262, 612            # crop box of panel b, PDF points (y down)
s = DPI / 72
subprocess.run(['pdftoppm', '-f', '4', '-l', '4', '-r', str(DPI), '-png', '-singlefile',
                '-x', str(int(X0*s)), '-y', str(int(Y0*s)), '-W', str(int((X1-X0)*s)), '-H', str(int((Y1-Y0)*s)),
                PDF, 'fig2b_crop'], check=True)
img = Image.open('fig2b_crop.png')

# axis calibration from tick-label centres (pdftotext -bbox), PDF points
DEC_Y = 24.33                                  # pt per decade on y (all panels)
y_at_CV1 = {0: 301.12, 1: 417.10, 2: 530.60}   # row -> y of CV = 1
count_cols = {0: (132.67, 21.05), 1: (276.20, 20.87), 2: (420.10, 20.93)}  # col -> (x at mean 0.1, pt/decade)
panels = {(r, c): count_cols[c] for r in range(3) for c in range(3)}
panels[(1, 2)] = ('quotient', 407.92, 23.77)   # Sanity: x at 1e-6
panels[(2, 2)] = ('quotient', 407.92, 23.84)   # scVI

fig = plt.figure(figsize=((X1-X0)/72*1.6, (Y1-Y0)/72*1.6))
ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off')
ax.imshow(img, extent=[X0, X1, Y1, Y0])
ax.set_xlim(X0, X1); ax.set_ylim(Y1, Y0)

for (r, c), cal in panels.items():
    yc1 = y_at_CV1[r]
    ypt = lambda cv: yc1 - DEC_Y * np.log10(cv)
    if cal[0] == 'quotient':                   # x in transcription quotients: slope only, position arbitrary
        _, x6, dx = cal
        xpt = lambda m: x6 + dx * np.log10(m / 1e-6)
        m = np.logspace(-6.3, -2, 50); cv = 0.02 * (m / 1e-3) ** -0.5
        style = dict(color='black', lw=1.4, ls='--')
        label = 'slope $-1/2$'
    else:                                      # x in UMI-scale counts: exact Poisson CV = 1/sqrt(mean)
        x01, dx = cal
        xpt = lambda m: x01 + dx * np.log10(m / 0.1)
        m = np.logspace(-2, 3, 50); cv = m ** -0.5
        style = dict(color='black', lw=1.4)
        label = r'Poisson: CV$=1/\sqrt{\rm mean}$'
    keep = (cv <= 10) & (cv >= 0.006)
    ax.plot(xpt(m[keep]), ypt(cv[keep]), **style)
    i = int(keep.sum() * 0.62); xs, ys = xpt(m[keep])[i], ypt(cv[keep])[i]
    ax.text(xs + 2, ys + 2.5, label, fontsize=6.5, va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.8))

fig.savefig('Fig2b_Poisson_slope.pdf', dpi=DPI)
fig.savefig('Fig2b_Poisson_slope.png', dpi=200)
