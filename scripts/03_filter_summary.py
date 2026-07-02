"""
DESeq2 sonuçlarını (results/deseq2_results.csv) anlamlılık eşiklerine göre
süzer, özet istatistik yazdırır ve volkan/MA grafiklerini üretir.

Eşik seçimi (padj < 0.05 ve |log2FoldChange| > 1): padj < 0.05, yanlış
pozitif oranını (Benjamini-Hochberg ile düzeltilmiş) kontrol eder;
log2FoldChange eşiği ise sadece istatistiksel olarak anlamlı değil, aynı
zamanda biyolojik olarak da göz ardı edilemeyecek büyüklükte (>2 kat)
bir değişimi olan genleri tutmak içindir — çok küçük ama istatistiksel
olarak "anlamlı" etkiler (ör. log2FC=0.1) büyük örneklem/derin
sekanslamada kolayca anlamlı çıkabilir ama biyolojik olarak ilgi çekici
olmayabilir.
"""

import matplotlib

matplotlib.use("Agg")  # headless ortam: görüntü sunucusu yok, doğrudan PNG'ye yaz

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
DESEQ_RESULTS_PATH = RESULTS_DIR / "deseq2_results.csv"
SIGNIFICANT_GENES_PATH = RESULTS_DIR / "significant_genes.csv"

PADJ_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1.0

# dataviz stilinden: log2FoldChange doğası gereği diverging (sıfır etrafında
# yön taşıyan) bir ölçü, bu yüzden mavi/kırmızı diverging çift kullanılıyor;
# anlamlı olmayan genler nötr gri.
COLOR_NS = "#c3c2b7"
COLOR_UP = "#e34948"
COLOR_DOWN = "#2a78d6"
COLOR_MUTED_LINE = "#898781"
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e1e0d9"


def load_results() -> pd.DataFrame:
    df = pd.read_csv(DESEQ_RESULTS_PATH, index_col=0)
    df.index.name = "gene_id"
    return df


def classify(df: pd.DataFrame) -> pd.Series:
    """Her gen için 'up' / 'down' / 'ns' etiketi üretir. padj NaN olan
    genler (DESeq2'nin independent filtering'i düşük count'lu genler
    için padj'ı NA bırakır) anlamlı sayılmaz, sessizce 'ns' olur."""
    is_sig = (df["padj"] < PADJ_THRESHOLD) & (df["log2FoldChange"].abs() > LOG2FC_THRESHOLD)
    is_sig = is_sig.fillna(False)
    direction = np.where(df["log2FoldChange"] > 0, "up", "down")
    return pd.Series(np.where(is_sig, direction, "ns"), index=df.index)


def summarize(df: pd.DataFrame, label: pd.Series) -> None:
    n_total = len(df)
    n_up = (label == "up").sum()
    n_down = (label == "down").sum()
    n_sig = n_up + n_down

    print(f"Toplam test edilen gen: {n_total}")
    print(
        f"Anlamlı gen (padj < {PADJ_THRESHOLD} ve |log2FC| > {LOG2FC_THRESHOLD}): {n_sig}"
    )
    print(f"  yukarı (log2FC > 0): {n_up}")
    print(f"  aşağı  (log2FC < 0): {n_down}")

    print("\nEn anlamlı 10 gen (padj'a göre):")
    top10 = df.sort_values("padj").head(10)
    print(top10[["baseMean", "log2FoldChange", "padj"]])


def write_significant_genes(df: pd.DataFrame, label: pd.Series) -> pd.DataFrame:
    sig_df = df[label != "ns"].copy()
    sig_df["direction"] = label[label != "ns"]
    sig_df = sig_df.sort_values("padj")
    sig_df.to_csv(SIGNIFICANT_GENES_PATH)
    print(f"\n[tamam] anlamlı genler yazıldı: {SIGNIFICANT_GENES_PATH} ({len(sig_df)} satır)")
    return sig_df


def _style_axes(ax) -> None:
    """Ortak eksen görünümü: nötr grid, resesif spine'lar — verinin
    kendisi öne çıksın diye çerçeve/gridline'lar sessiz tutuluyor."""
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=INK_SECONDARY, labelsize=9)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(COLOR_MUTED_LINE)
    ax.grid(True, color=GRID, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)


def _scatter_by_label(ax, x: pd.Series, y: pd.Series, label: pd.Series) -> None:
    """ns -> down -> up sırasıyla çiziyoruz ki anlamlı noktalar
    (down/up) görsel olarak arka plandaki gri bulutun üzerinde kalsın."""
    for key, color, legend_label in (
        ("ns", COLOR_NS, "anlamlı değil"),
        ("down", COLOR_DOWN, "aşağı"),
        ("up", COLOR_UP, "yukarı"),
    ):
        mask = label == key
        ax.scatter(
            x[mask],
            y[mask],
            s=10 if key == "ns" else 14,
            color=color,
            alpha=0.5 if key == "ns" else 0.85,
            linewidths=0,
            label=f"{legend_label} (n={mask.sum()})",
        )


def plot_volcano(df: pd.DataFrame, label: pd.Series, dest: Path) -> None:
    """x = log2FoldChange, y = -log10(padj)."""
    plot_df = df.dropna(subset=["padj", "log2FoldChange"]).copy()
    plot_label = label.loc[plot_df.index]

    # padj tam olarak 0 olursa -log10(0) sonsuz olur; pratikte DESeq2
    # bunu üretmez ama garanti altına almak için en küçük pozitif padj
    # değerine clip ediyoruz.
    positive_padj = plot_df.loc[plot_df["padj"] > 0, "padj"]
    padj_floor = positive_padj.min() if len(positive_padj) else 1e-300
    plot_df["neglog10_padj"] = -np.log10(plot_df["padj"].clip(lower=padj_floor))

    fig, ax = plt.subplots(figsize=(7, 6), facecolor=SURFACE)
    _style_axes(ax)
    _scatter_by_label(ax, plot_df["log2FoldChange"], plot_df["neglog10_padj"], plot_label)

    ax.axvline(LOG2FC_THRESHOLD, color=COLOR_MUTED_LINE, linestyle="--", linewidth=1)
    ax.axvline(-LOG2FC_THRESHOLD, color=COLOR_MUTED_LINE, linestyle="--", linewidth=1)
    ax.axhline(-np.log10(PADJ_THRESHOLD), color=COLOR_MUTED_LINE, linestyle="--", linewidth=1)

    ax.set_xlabel("log2 fold change (dexamethasone / untreated)", color=INK_SECONDARY)
    ax.set_ylabel("-log10(padj)", color=INK_SECONDARY)
    ax.set_title("Volkan grafiği: dexamethasone vs untreated", color=INK_PRIMARY, fontsize=13, loc="left")
    ax.legend(frameon=False, loc="upper left", fontsize=9, labelcolor=INK_SECONDARY)

    fig.tight_layout()
    fig.savefig(dest, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"[tamam] volkan grafiği yazıldı: {dest}")


def plot_ma(df: pd.DataFrame, label: pd.Series, dest: Path) -> None:
    """x = log2(baseMean + 1), y = log2FoldChange. baseMean sıfır
    olabildiği için log2(0) tanımsızlığını önlemek üzere +1 ekleniyor."""
    plot_df = df.dropna(subset=["log2FoldChange", "baseMean"]).copy()
    plot_df["log2_basemean"] = np.log2(plot_df["baseMean"] + 1)
    plot_label = label.loc[plot_df.index]

    fig, ax = plt.subplots(figsize=(7, 6), facecolor=SURFACE)
    _style_axes(ax)
    _scatter_by_label(ax, plot_df["log2_basemean"], plot_df["log2FoldChange"], plot_label)

    ax.axhline(0, color=COLOR_MUTED_LINE, linestyle="-", linewidth=1)

    ax.set_xlabel("log2(baseMean + 1)", color=INK_SECONDARY)
    ax.set_ylabel("log2 fold change (dexamethasone / untreated)", color=INK_SECONDARY)
    ax.set_title("MA grafiği: dexamethasone vs untreated", color=INK_PRIMARY, fontsize=13, loc="left")
    ax.legend(frameon=False, loc="upper right", fontsize=9, labelcolor=INK_SECONDARY)

    fig.tight_layout()
    fig.savefig(dest, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"[tamam] MA grafiği yazıldı: {dest}")


def main() -> None:
    df = load_results()
    label = classify(df)

    summarize(df, label)
    write_significant_genes(df, label)

    plot_volcano(df, label, RESULTS_DIR / "volcano_plot.png")
    plot_ma(df, label, RESULTS_DIR / "ma_plot.png")


if __name__ == "__main__":
    main()
