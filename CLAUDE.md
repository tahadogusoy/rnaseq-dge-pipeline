# CLAUDE.md

Bu dosya, bu proje üzerinde çalışırken Claude Code'un uyması gereken bağlamı
ve kuralları tanımlar.

## Proje

Bu proje, HSC (hematopoietik kök hücre) ve MSC (mezenkimal kök hücre)
örnekleri arasında bulk RNA-seq verisi kullanarak diferansiyel gen ifadesi
(differential gene expression) analizi yapan bir pipeline'dır.

## Amaç

- Ham count verisinden başlayarak HSC vs MSC karşılaştırmasında anlamlı
  şekilde farklı ifade edilen genleri bulmak
- Sonuçları DESeq2 ile istatistiksel olarak değerlendirmek
- clusterProfiler ile GO/KEGG zenginleştirme analizi yaparak biyolojik
  yorum katmanı eklemek

## Ortam

- Conda ortamı: `dge` (bkz. `environment.yml`)
- Python 3.11 (pandas, numpy) veri hazırlama ve yardımcı işler için
- R (r-base, DESeq2, clusterProfiler) istatistiksel diferansiyel ifade ve
  zenginleştirme analizi için
- Ortamı aktive etmeden hiçbir script çalıştırılmamalı:
  `conda activate dge`

## scripts/ Mantığı

Betikler `scripts/` altında numaralandırılarak sırayla çalışır ve her biri
bir öncekinin çıktısını girdi olarak alır:

- `01_...` : ham veriyi (`data/raw/`) okuma, temizleme, `data/processed/`'a yazma
- `02_...` : DESeq2 ile diferansiyel ifade analizi
- `03_...` : sonuçların filtrelenmesi/özetlenmesi
- `04_...` : clusterProfiler ile GO/KEGG zenginleştirme analizi

Her betik bağımsız çalıştırılabilir olmalı ama sıralamaya (01 → 02 → 03 → 04)
uyulduğu varsayılır. Yeni bir adım eklenirse bir sonraki numarayla eklenmeli,
mevcut sıralama bozulmamalı. Ara çıktılar `data/processed/` içine, nihai
analiz çıktıları (grafik, tablo, rapor) `results/` içine yazılır.

## Kod Yazma Kuralı

Kod yazarken (özellikle istatistiksel/biyoinformatik adımlarda: filtreleme
eşikleri, normalizasyon yöntemi, test seçimi, parametre seçimleri) neden o
seçimin yapıldığı açıklanmalı. Sadece "ne yapıldığı" değil, "neden öyle
yapıldığı" da yorum satırlarında veya commit/PR açıklamalarında belirtilmeli.
Bu, biyolojik/istatistiksel varsayımların sonradan sorgulanabilir ve
denetlenebilir olmasını sağlar.
