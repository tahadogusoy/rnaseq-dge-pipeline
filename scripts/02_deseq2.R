# DESeq2 ile diferansiyel gen ifadesi analizi: dexamethasone vs untreated.
#
# Neden bu karşılaştırma: data/processed/sample_metadata.csv içinde 4
# treatment grubu var (Untreated, Dexamethasone, Albuterol,
# Albuterol_Dexamethasone). Albuterol kombinasyonlarını dahil etmek
# "condition" etkisini bir başka değişkenle (albuterol) karıştırır;
# temiz bir "dexamethasone etkisi" tahmini için sadece Untreated
# (referans) ve Dexamethasone örnekleri tutuluyor.
#
# Neden design = ~ cellline + condition: örnekler 4 hücre hattından
# (cell line) tekrarlı alınmış (her hücre hattı hem Untreated hem
# Dexamethasone içeriyor). Hücre hatları arasında bazal ifade farkları
# olması beklenir; bunu kovaryat olarak modele eklemezsek, hücre
# hattından kaynaklanan varyans hataya (residual) karışır ve
# dexamethasone etkisinin gücü (istatistiksel power) azalır.

suppressPackageStartupMessages(library(DESeq2))

raw_dir <- file.path("data", "raw")
processed_dir <- file.path("data", "processed")
results_dir <- file.path("results")

# Accession'ı sabitlemek yerine data/raw altında dosya adına göre arıyoruz,
# çünkü 01_download.py --gse parametrik: hangi GSE indirildiyse o
# kullanılsın istiyoruz.
count_files <- list.files(
  raw_dir,
  pattern = "_raw_counts_GRCh38\\.p13_NCBI\\.tsv\\.gz$",
  full.names = TRUE
)
if (length(count_files) == 0) {
  stop(
    "data/raw altında ham count dosyası bulunamadı. ",
    "Önce 'conda run -n dge python scripts/01_download.py' çalıştırılmalı."
  )
}
if (length(count_files) > 1) {
  stop(
    "data/raw altında birden fazla ham count dosyası var, hangisi ",
    "kullanılacak belirsiz: ", paste(basename(count_files), collapse = ", ")
  )
}
count_file <- count_files[1]
metadata_file <- file.path(processed_dir, "sample_metadata.csv")
if (!file.exists(metadata_file)) {
  stop(
    "sample_metadata.csv bulunamadı (", metadata_file, "). ",
    "Önce 'conda run -n dge python scripts/01_download.py' çalıştırılmalı."
  )
}

cat("Ham count dosyası:", count_file, "\n")
counts <- read.delim(gzfile(count_file), row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
# 01_download.py bu matrisin tam sayı olduğunu zaten doğruladı; burada
# R tarafında da tipi integer'a çeviriyoruz çünkü DESeq2 negatif binom
# modeli integer count bekler.
storage.mode(counts) <- "integer"

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
if (!all(metadata$sample_id %in% colnames(counts))) {
  missing_ids <- setdiff(metadata$sample_id, colnames(counts))
  stop(
    "sample_metadata.csv içindeki bazı sample_id'ler count matrisinde yok: ",
    paste(missing_ids, collapse = ", ")
  )
}

cat("\n=== Filtrelemeden ÖNCE treatment etiketleri ve grup boyutları ===\n")
print(table(metadata$treatment))

# Beklenen etiketleri sessizce yok saymak yerine açıkça doğruluyoruz:
# GEO'daki characteristics alanları serbest metin olduğu için farklı bir
# GSE ile çalıştırıldığında (--gse) bu etiketler hiç bulunmayabilir.
required_levels <- c("Untreated", "Dexamethasone")
missing_levels <- setdiff(required_levels, unique(metadata$treatment))
if (length(missing_levels) > 0) {
  stop(
    "Beklenen treatment etiketleri metadata'da bulunamadı: ",
    paste(missing_levels, collapse = ", "),
    ". Mevcut etiketler: ", paste(unique(metadata$treatment), collapse = ", ")
  )
}

metadata_sub <- metadata[metadata$treatment %in% required_levels, ]
# required_levels sadece "Untreated" ve "Dexamethasone" içerdiği için bu
# filtre albuterol içeren tüm grupları (Albuterol, Albuterol_Dexamethasone)
# otomatik olarak dışarıda bırakıyor; yine de varsayımı sessizce
# güvenmek yerine açıkça doğruluyoruz.
if (any(grepl("Albuterol", metadata_sub$treatment, ignore.case = TRUE))) {
  stop("Filtreleme sonrası albuterol içeren örnekler kaldı; beklenmeyen durum.")
}

cat("\n=== Filtrelemeden SONRA kullanılacak örnekler ===\n")
print(table(metadata_sub$treatment))

counts_sub <- counts[, metadata_sub$sample_id, drop = FALSE]

col_data <- metadata_sub
rownames(col_data) <- col_data$sample_id
col_data <- col_data[colnames(counts_sub), , drop = FALSE]  # sıralamayı garanti et
# "Untreated" referans (ilk) seviye olacak şekilde elle sıralıyoruz;
# aksi halde factor() alfabetik sıralar ("Dexamethasone" < "Untreated")
# ve referans yanlışlıkla Dexamethasone olurdu.
col_data$condition <- factor(col_data$treatment, levels = required_levels)
col_data$cellline <- factor(col_data$cellline)

dds <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData = col_data,
  design = ~ cellline + condition
)

# DESeq2 vignette'inin önerdiği standart bir ön-filtre: toplam count'u
# çok düşük (<10) olan genler modele neredeyse hiç bilgi katmaz, sadece
# hesaplama yükünü ve çoklu test (BH/padj) düzeltmesindeki güç kaybını
# artırır. Bu, results() içindeki otomatik bağımsız filtrelemeden
# farklı, modelleme öncesi kaba bir ön temizlik.
keep <- rowSums(counts(dds)) >= 10
cat(sprintf(
  "\nÖn-filtre: %d / %d gen tutuldu (toplam count >= 10)\n",
  sum(keep), length(keep)
))
dds <- dds[keep, ]

dds <- DESeq(dds)

# contrast'ı açıkça belirtiyoruz: design'ın son terimi 'condition' olsa
# da results()'ın varsayılan yönü hangi seviyenin referans olduğuna
# bağlıdır; burada niyeti (dexamethasone / untreated) koddan okunur
# kılmak için elle yazıyoruz.
res <- results(dds, contrast = c("condition", "Dexamethasone", "Untreated"))
res_ordered <- res[order(res$padj), ]

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(results_dir, "deseq2_results.csv")
write.csv(as.data.frame(res_ordered), file = out_file)
cat("\n[tamam] sonuçlar yazıldı:", out_file, "\n")

cat("\n=== results() summary (dexamethasone vs untreated) ===\n")
summary(res)

cat("\n=== En anlamlı 10 gen (padj'a göre sıralı) ===\n")
print(head(as.data.frame(res_ordered), 10))
