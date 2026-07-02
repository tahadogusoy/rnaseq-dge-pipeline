# HSC vs MSC Bulk RNA-seq Diferansiyel İfade Analizi

Bu proje, hematopoietik kök hücreler (HSC) ve mezenkimal kök hücreler (MSC)
arasındaki gen ifade farklılıklarını bulk RNA-seq verisi üzerinden analiz
etmeyi amaçlar.

## Amaç

- Ham okuma sayılarından (raw counts) diferansiyel ifade analizi yapmak
- HSC ve MSC grupları arasında anlamlı farklılık gösteren genleri belirlemek
- Sonuçları fonksiyonel zenginleştirme (GO/KEGG) analizi ile yorumlamak

## Klasör Yapısı

```
data/raw/         # Ham veri (git'e girmez)
data/processed/   # İşlenmiş/temizlenmiş veri (git'e girmez)
scripts/          # Sıralı çalışan analiz betikleri (01 -> 02 -> 03 -> 04)
results/          # Analiz çıktıları, grafikler, tablolar (git'e girmez)
notebooks/         # Keşifsel analiz not defterleri
```

## Ortam Kurulumu

```bash
conda env create -f environment.yml
conda activate dge
```

## Pipeline

Analiz betikleri `scripts/` altında numaralandırılmış sırayla çalıştırılır,
her biri bir öncekinin çıktısını girdi olarak alır:

- `01_download.py`: GEO'dan ham count matrisini ve örnek metadata'sını indirir
  (`data/raw/`), örnek metadata'sını temizleyip `data/processed/` altına yazar
- `02_deseq2.R`: DESeq2 ile diferansiyel gen ifadesi analizi yapar
- `03_filter_summary.py`: sonuçları anlamlılık eşiklerine göre süzer, volkan
  ve MA grafiklerini üretir
- `04_enrichment.R`: clusterProfiler ile GO (Biological Process) ve KEGG
  yolak zenginleştirme analizi yapar

Tüm çıktılar (tablolar, grafikler) `results/` altına yazılır.

## Durum

Uçtan uca pipeline (01 → 04) çalışır durumda.
