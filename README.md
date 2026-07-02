# Bulk RNA-seq Diferansiyel Gen İfadesi ve Zenginleştirme Analizi Pipeline'ı

GEO'dan indirilen bulk RNA-seq ham count verisi üzerinde DESeq2 ile
diferansiyel gen ifadesi analizi yapan ve sonuçları clusterProfiler ile
GO/KEGG zenginleştirme analiziyle yorumlayan, uçtan uca tekrar
üretilebilir bir pipeline.

## Veri

Depodaki örnek çalıştırma **GSE52778** (insan havayolu düz kas hücreleri,
dexamethasone tedavisi) verisiyle yapılmıştır. Bu veri seti belirli bir
biyolojik hipotez için değil, pipeline'ın baştan sona doğru çalıştığını
göstermek için seçilmiştir: NCBI'nin GEO üzerinde her seri için ürettiği
standart ham count matrisine sahiptir ve iyi bilinen, temiz bir deney
tasarımı (tedavi vs kontrol) içerir.

Pipeline parametriktir; `scripts/01_download.py` bir `--gse` argümanı
alır ve herhangi bir GEO serisine yönlendirilebilir (script, o seri için
NCBI'nin ürettiği ham count matrisini otomatik olarak arar).

## Pipeline Adımları

| Adım | Script | Ne yapar |
|---|---|---|
| 01 | `scripts/01_download.py` | GEO'dan NCBI'nin ürettiği ham count matrisini ve örnek (sample) metadata'sını indirir, `data/processed/sample_metadata.csv`'yi üretir |
| 02 | `scripts/02_deseq2.R` | DESeq2 ile diferansiyel gen ifadesi analizini çalıştırır, sonuçları `results/deseq2_results.csv`'ye yazar |
| 03 | `scripts/03_filter_summary.py` | anlamlılık eşiklerine göre genleri filtreler, volkan ve MA grafiklerini üretir |
| 04 | `scripts/04_enrichment.R` | clusterProfiler ile GO (Biological Process) ve KEGG yolak zenginleştirme analizini yapar |

Her adım bir öncekinin çıktısını girdi olarak alır; ara ve nihai çıktılar
`results/` altına yazılır.

## Kullanılan Araçlar

- Python 3.11 (pandas, numpy, requests)
- R (DESeq2, clusterProfiler, org.Hs.eg.db)
- conda (`environment.yml` ile ortam yönetimi)

## Kurulum ve Çalıştırma

```bash
conda env create -f environment.yml
conda activate dge

conda run -n dge python scripts/01_download.py
conda run -n dge Rscript scripts/02_deseq2.R
conda run -n dge python scripts/03_filter_summary.py
conda run -n dge Rscript scripts/04_enrichment.R
```

## Örnek Sonuçlar (GSE52778: dexamethasone vs untreated)

- Ön-filtreden sonra **21.833 gen** test edildi
- Anlamlılık eşiği (padj < 0.05 ve |log2FC| > 1) ile **1.157 anlamlı gen**:
  584 yukarı, 573 aşağı regüle
- GO (Biological Process) zenginleştirmesinde **588 anlamlı terim**
  (qvalue < 0.05)
- KEGG yolak zenginleştirmesinde **17 anlamlı yolak** (qvalue < 0.05)

## Not

Bu projedeki metodolojik kararlar (analiz tasarımı, veri/GEO serisi
seçimi, kovaryat ve kontrast tanımları, filtreleme eşikleri) tarafımca
alınmıştır. Kodlama aşamasında yapay zeka destekli bir asistandan
yararlanılmıştır.
