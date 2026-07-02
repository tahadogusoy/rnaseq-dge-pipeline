"""
GEO'dan bir çalışmanın NCBI ham (raw) count matrisini ve örnek (sample)
metadata'sını indirir, data/raw/ altına yazar ve indirilen veriyi
DESeq2 girdisine uygunluk açısından denetler.

Neden bu iki dosya:
- NCBI'nin "rnaseq_counts" formatı (GSE..._raw_counts_GRCh38.p13_NCBI.tsv.gz),
  GEO'ya yüklenen ham FASTQ'lardan NCBI tarafından yeniden, tek bir ortak
  pipeline (STAR + featureCounts benzeri) ile üretilmiş HAM count matrisidir.
  Yazarın kendi işlediği "processed" dosyalar genelde zaten normalize
  edilmiş olabiliyor (TPM/FPKM/DESeq2 normalized) ve DESeq2'ye ham count
  olmadan verilemez; bu yüzden yazarın dosyası yerine NCBI'nin kendi
  yeniden işlediği ham count dosyası tercih edilir.
- series_matrix.txt.gz dosyası ise örneklerin (GSM) başlık ve
  characteristics (grup, doku, hasta/kontrol vb.) bilgisini standart bir
  formatta içerir; bu bilgi olmadan count matrisindeki sütunları hangi
  deney grubuna ait olduğunu bilemeyiz.
"""

import argparse
import csv
import gzip
import sys
import urllib.error
from collections import Counter
from pathlib import Path

import requests

RAW_COUNTS_URL_TEMPLATE = (
    "https://www.ncbi.nlm.nih.gov/geo/download/"
    "?type=rnaseq_counts&acc={gse}&format=file"
    "&file={gse}_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)
SERIES_MATRIX_URL_TEMPLATE = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_folder}nnn/{gse}/matrix/"
    "{gse}_series_matrix.txt.gz"
)

DATA_RAW_DIR = Path(__file__).resolve().parent.parent / "data" / "raw"
DATA_PROCESSED_DIR = Path(__file__).resolve().parent.parent / "data" / "processed"


def gse_folder(gse: str) -> str:
    """GSE52778 -> 'GSE52' (GEO'nun FTP dizin şeması binlik değil, son 3
    hane atılarak oluşturuluyor: GSE{n}nnn)."""
    digits = gse.removeprefix("GSE")
    return "GSE" + digits[:-3]


def download_file(url: str, dest: Path, description: str) -> Path:
    """Tek bir dosyayı indirir. 404/erişim hatalarını sessizce yutmaz,
    çünkü ham count dosyası olmadan pipeline'ın geri kalanı anlamsız
    sonuç üretir ama fark edilmeden devam edebilir."""
    print(f"[indiriliyor] {description}: {url}")
    response = requests.get(url, timeout=60)
    if response.status_code == 404:
        raise FileNotFoundError(
            f"HATA: {description} bulunamadı (404). "
            f"Bu GSE için NCBI'nin ürettiği ham count matrisi yok olabilir. "
            f"URL: {url}"
        )
    response.raise_for_status()
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(response.content)
    print(f"[tamam] {dest.name} ({dest.stat().st_size:,} bayt)")
    return dest


def inspect_count_matrix(path: Path) -> None:
    """Count matrisinin boyutunu, ilk satırlarını yazdırır ve DESeq2'nin
    gerektirdiği gibi ham (integer) count olup olmadığını kontrol eder.
    DESeq2 negatif binom modeli tam sayı count varsayar; ondalıklı
    (normalize/TPM) değerler verilirse model varsayımı bozulur ve
    sonuçlar istatistiksel olarak geçersiz olur — bu yüzden bu kontrol
    atlanamaz."""
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in f]

    n_genes = len(rows)
    n_samples = len(header) - 1  # ilk kolon gen ID'si
    print(f"\n=== Count matrisi: {path.name} ===")
    print(f"Satır (gen) sayısı : {n_genes}")
    print(f"Sütun (örnek) sayısı: {n_samples}")
    print("\nİlk 5 satır:")
    print("\t".join(header))
    for row in rows[:5]:
        print("\t".join(row))

    # Tüm değerleri değil, örnekleme yaparak kontrol ediyoruz: matris
    # büyük olabilir, tam tarama gereksiz maliyetli ve amaç zaten sadece
    # "ham count mu değil mi" sorusuna hızlı bir cevap vermek.
    sample_values = [v for row in rows[:200] for v in row[1:]]
    non_integer = [v for v in sample_values if not _is_integer_string(v)]

    print("\n=== Ham count kontrolü ===")
    if non_integer:
        print(
            f"UYARI: Değerler TAM SAYI DEĞİL (örnek: {non_integer[:5]}). "
            "Bu muhtemelen normalize edilmiş (TPM/FPKM/DESeq2-normalized) "
            "bir matris ve DESeq2'ye HAM count olarak verilemez."
        )
    else:
        print(
            "OK: örneklenen değerlerin tümü tam sayı — ham count matrisiyle "
            "tutarlı (DESeq2 girdisi için uygun görünüyor)."
        )


def _is_integer_string(value: str) -> bool:
    try:
        float_val = float(value)
    except ValueError:
        return False
    return float_val.is_integer()


def parse_series_matrix(path: Path) -> list[dict]:
    """series_matrix.txt.gz içinden her örnek için sample_id (GSM), title
    ve characteristics alanlarını çıkarır.

    Characteristics satırlarının sırası (treatment / tissue / cell line /
    ...) GEO'da yazarın yükleme sırasına bağlıdır ve GSE'den GSE'ye
    değişebilir; bu yüzden sütun indeksine değil, her hücredeki
    'key: value' içindeki key'e göre eşleştiriyoruz — sıra değişse bile
    doğru çalışır."""
    titles: list[str] = []
    sample_ids: list[str] = []
    characteristics_rows: list[list[str]] = []

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                titles = _split_series_matrix_line(line)
            elif line.startswith("!Sample_geo_accession"):
                sample_ids = _split_series_matrix_line(line)
            elif line.startswith("!Sample_characteristics_ch1"):
                characteristics_rows.append(_split_series_matrix_line(line))

    if not titles or not sample_ids:
        raise ValueError(
            "series matrix beklenen formatta değil: !Sample_title veya "
            "!Sample_geo_accession satırı bulunamadı."
        )

    samples = []
    for i, (sample_id, title) in enumerate(zip(sample_ids, titles)):
        char_dict = {}
        for row in characteristics_rows:
            if i >= len(row):
                continue
            key, _, value = row[i].partition(": ")
            char_dict[key.strip().lower()] = value.strip()
        samples.append(
            {
                "sample_id": sample_id,
                "title": title,
                "cellline": char_dict.get("cell line", ""),
                "treatment": char_dict.get("treatment", ""),
            }
        )
    return samples


def inspect_series_matrix(samples: list[dict]) -> None:
    print(f"\n=== Örnek metadata ({len(samples)} örnek) ===")
    for s in samples:
        print(
            f"- {s['sample_id']} ({s['title']}): "
            f"cellline={s['cellline']}, treatment={s['treatment']}"
        )

    # 02_deseq2.R hangi grupları tutup hangilerini hariç tutacağına karar
    # vermeden önce bu dağılımı görmek gerekiyor (örn. albuterol
    # kombinasyonlarının kaç örnek olduğu).
    print("\n=== Treatment grup dağılımı ===")
    counts = Counter(s["treatment"] for s in samples)
    for treatment, n in sorted(counts.items()):
        print(f"{treatment}: {n} örnek")


def write_sample_metadata(samples: list[dict], dest: Path) -> None:
    """Örnek metadata'sını data/processed/sample_metadata.csv olarak yazar.
    sample_id değerleri, count matrisindeki GSM sütun adlarıyla birebir
    aynı olacak şekilde series matrix'ten alınıyor; 02_deseq2.R bu
    dosyayı count matrisi sütunlarıyla eşleştirip colData oluşturmak
    için kullanacak."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["sample_id", "title", "cellline", "treatment"]
        )
        writer.writeheader()
        writer.writerows(samples)
    print(f"\n[tamam] örnek metadata yazıldı: {dest} ({len(samples)} satır)")


def _split_series_matrix_line(line: str) -> list[str]:
    """'!Sample_title\t"a"\t"b"\n' -> ['a', 'b']"""
    parts = line.rstrip("\n").split("\t")[1:]
    return [p.strip('"') for p in parts]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="GEO'dan ham count matrisi ve örnek metadata indirir."
    )
    parser.add_argument(
        "--gse",
        default="GSE52778",
        help="GEO series accession (örn. GSE52778). Varsayılan: GSE52778.",
    )
    args = parser.parse_args()
    gse = args.gse

    raw_counts_url = RAW_COUNTS_URL_TEMPLATE.format(gse=gse)
    series_matrix_url = SERIES_MATRIX_URL_TEMPLATE.format(
        gse_folder=gse_folder(gse), gse=gse
    )

    raw_counts_path = DATA_RAW_DIR / f"{gse}_raw_counts_GRCh38.p13_NCBI.tsv.gz"
    series_matrix_path = DATA_RAW_DIR / f"{gse}_series_matrix.txt.gz"

    try:
        download_file(raw_counts_url, raw_counts_path, "NCBI ham count matrisi")
        download_file(series_matrix_url, series_matrix_path, "series matrix (örnek metadata)")
    except (FileNotFoundError, requests.HTTPError, requests.ConnectionError) as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)

    print("\n=== İndirilen dosyalar ===")
    for p in (raw_counts_path, series_matrix_path):
        print(f"{p.name}: {p.stat().st_size:,} bayt")

    inspect_count_matrix(raw_counts_path)

    samples = parse_series_matrix(series_matrix_path)
    inspect_series_matrix(samples)
    write_sample_metadata(samples, DATA_PROCESSED_DIR / "sample_metadata.csv")


if __name__ == "__main__":
    main()
