# GO (Biological Process) ve KEGG yolak zenginleştirme analizi.
#
# Neden ENTREZ ID: NCBI'nin ürettiği ham count matrisindeki GeneID sütunu
# NCBI Gene (Entrez) numaralarıdır (ör. "783", "8404"); org.Hs.eg.db ve
# clusterProfiler'ın enrichGO/enrichKEGG fonksiyonları bu ID uzayını
# doğrudan (keyType="ENTREZID") destekliyor, bu yüzden ayrıca bir ID
# dönüşümüne (ör. Symbol/Ensembl) gerek yok — ama yine de ID'lerin
# gerçekten sayısal Entrez formatında olduğunu doğruluyoruz, çünkü
# CSV'den okunan bir sütunun sessizce bozuk (boşluk, tırnak, NA) gelmesi
# ihtimali var ve bu durumu fark etmeden enrichGO'ya vermek, genlerin
# hiç eşleşmemesine ve "0 gen tutuldu"nun sessizce geçmesine yol açabilir.
#
# Neden universe = tüm test edilen genler: clusterProfiler'a universe
# verilmezse varsayılan arka plan tüm insan genomudur. Bizim deneyimiz
# ise ön-filtreden (rowSums >= 10) geçen ~22 bin genle sınırlı; gerçek
# arka plan olarak tüm genomu kullanmak, bu dokuda zaten yüksek ifade
# edilen / test edilebilir olan yolakları yapay olarak "zenginleşmiş"
# gösterir (seçim yanlılığı). Bu yüzden universe'i deseq2_results.csv'deki
# TÜM test edilen genler olarak açıkça veriyoruz.

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

# Bazı sistemlerde ggsave/png() ekran sunucusu olmadan varsayılan
# bitmap tipiyle sorun çıkarabiliyor; cairo headless ortamlarda da
# güvenilir çalışıyor.
options(bitmapType = "cairo")

results_dir <- file.path("results")
sig_genes_path <- file.path(results_dir, "significant_genes.csv")
all_results_path <- file.path(results_dir, "deseq2_results.csv")

if (!file.exists(sig_genes_path) || !file.exists(all_results_path)) {
  stop(
    "results/significant_genes.csv veya results/deseq2_results.csv yok. ",
    "Önce 'conda run -n dge python scripts/03_filter_summary.py' çalıştırılmalı."
  )
}

sig_genes <- read.csv(sig_genes_path, stringsAsFactors = FALSE)
all_results <- read.csv(all_results_path, stringsAsFactors = FALSE)
# deseq2_results.csv'nin ilk sütunu R'ın write.csv'sinden kalma isimsiz
# satır adı (gene ID) sütunudur.
colnames(all_results)[1] <- "gene_id"

clean_entrez_ids <- function(raw_ids, label) {
  ids <- trimws(as.character(raw_ids))
  valid <- grepl("^[0-9]+$", ids)
  n_dropped <- sum(!valid)
  if (n_dropped > 0) {
    cat(sprintf(
      "UYARI: %s içinde %d ID geçerli bir Entrez ID formatında değil ve atıldı. Örnek: %s\n",
      label, n_dropped, paste(head(ids[!valid], 5), collapse = ", ")
    ))
  }
  unique(ids[valid])
}

sig_entrez <- clean_entrez_ids(sig_genes$gene_id, "significant_genes.csv")
universe_entrez <- clean_entrez_ids(all_results$gene_id, "deseq2_results.csv")

cat(sprintf(
  "\nGeçerli Entrez ID sayısı: anlamlı genler = %d, universe (tüm test edilenler) = %d\n",
  length(sig_entrez), length(universe_entrez)
))

# "Sessizce geçme, uyar" kuralı: enrichGO/enrichKEGG'e anlamlı sayıda
# gen vermeden çağrı yapmak istatistiksel olarak neredeyse hiçbir zaman
# anlamlı bir zenginleştirme üretmez ve bu durumu sonuç tablosuna
# bakmadan fark etmek zor olur.
MIN_SIG_GENES <- 10
if (length(sig_entrez) < MIN_SIG_GENES) {
  cat(sprintf(
    "UYARI: yalnızca %d anlamlı gen var (eşik: %d). Zenginleştirme sonuçları güvenilir olmayabilir.\n",
    length(sig_entrez), MIN_SIG_GENES
  ))
}
if (!all(sig_entrez %in% universe_entrez)) {
  n_outside <- sum(!(sig_entrez %in% universe_entrez))
  cat(sprintf(
    "UYARI: anlamlı genlerden %d tanesi universe listesinde yok (beklenmiyordu).\n",
    n_outside
  ))
}

report_enrichment <- function(result, label, n_top = 10) {
  df <- as.data.frame(result)
  if (nrow(df) == 0) {
    cat(sprintf("\nUYARI: %s için anlamlı (qvalue < 0.05) hiçbir terim bulunamadı.\n", label))
    return(df)
  }
  cat(sprintf("\n=== En anlamlı %d %s (qvalue'ya göre) ===\n", min(n_top, nrow(df)), label))
  print(head(df[order(df$qvalue), c("ID", "Description", "GeneRatio", "qvalue", "Count")], n_top))
  df
}

save_dotplot <- function(result, dest, title, n_top = 10) {
  df <- as.data.frame(result)
  if (nrow(df) == 0) {
    cat(sprintf("UYARI: %s boş olduğu için dotplot atlandı: %s\n", title, dest))
    return(invisible(NULL))
  }
  p <- tryCatch(
    dotplot(result, showCategory = n_top) + ggtitle(title),
    error = function(e) {
      cat(sprintf("UYARI: %s için dotplot üretilemedi: %s\n", title, conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(p)) {
    ggsave(dest, plot = p, width = 8, height = 6, dpi = 150, bg = "white")
    cat(sprintf("[tamam] dotplot yazıldı: %s\n", dest))
  }
}

## 1) GO Biological Process zenginleştirme -----------------------------
cat("\n--- GO (BP) zenginleştirme çalıştırılıyor ---\n")
ego <- enrichGO(
  gene = sig_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE  # gen sembollerini de CSV'ye ekler, Entrez ID'den daha okunur
)

go_out <- file.path(results_dir, "go_enrichment.csv")
write.csv(as.data.frame(ego), file = go_out, row.names = FALSE)
cat(sprintf("[tamam] GO sonuçları yazıldı: %s (%d terim)\n", go_out, nrow(as.data.frame(ego))))
report_enrichment(ego, "GO terimi")
save_dotplot(ego, file.path(results_dir, "go_dotplot.png"), "GO Biological Process zenginleştirme")

## 2) KEGG yolak zenginleştirme -----------------------------------------
# enrichKEGG, KEGG'in kendi REST API'sine canlı bir ağ isteği atar
# (organism="hsa" -> "Homo sapiens" KEGG kodu); internet yoksa veya
# KEGG servisi yanıt vermezse hata fırlatır. Bu hatayı sessizce yutmak
# yerine yakalayıp açıkça raporluyoruz, ama GO sonuçlarını (yukarıda
# zaten yazıldı) etkilememesi için script'in tamamını durdurmuyoruz.
cat("\n--- KEGG zenginleştirme çalıştırılıyor (KEGG REST API'sine bağlanıyor) ---\n")
ekegg <- tryCatch(
  enrichKEGG(
    gene = sig_entrez,
    universe = universe_entrez,
    organism = "hsa",
    keyType = "kegg",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  ),
  error = function(e) {
    cat(sprintf(
      "HATA: KEGG zenginleştirme başarısız oldu (muhtemelen ağ/KEGG API erişim sorunu): %s\n",
      conditionMessage(e)
    ))
    NULL
  }
)

if (is.null(ekegg)) {
  cat("UYARI: KEGG sonuçları üretilemedi, kegg_enrichment.csv yazılmadı.\n")
} else {
  ekegg_readable <- tryCatch(
    setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
    error = function(e) ekegg  # gen sembol dönüşümü başarısız olursa Entrez ID'lerle devam et
  )
  kegg_out <- file.path(results_dir, "kegg_enrichment.csv")
  write.csv(as.data.frame(ekegg_readable), file = kegg_out, row.names = FALSE)
  cat(sprintf(
    "[tamam] KEGG sonuçları yazıldı: %s (%d yolak)\n",
    kegg_out, nrow(as.data.frame(ekegg_readable))
  ))
  report_enrichment(ekegg_readable, "KEGG yolağı")
  save_dotplot(ekegg_readable, file.path(results_dir, "kegg_dotplot.png"), "KEGG yolak zenginleştirme")
}
