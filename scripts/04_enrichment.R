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

universe_entrez <- clean_entrez_ids(all_results$gene_id, "deseq2_results.csv")

# Neden yön (up/down) ayrı çalıştırılıyor: significant_genes.csv hem
# yukarı hem aşağı regüle edilen genleri "direction" sütunuyla ayırt
# edilebilir biçimde içeriyor (bkz. 03_filter_summary.py). Bunları tek
# bir zenginleştirme listesinde birleştirmek, dexamethasone tarafından
# İNDÜKLENEN ve BASKILANAN genlerin farklı biyolojik programları (ör.
# glukokortikoid reseptörünün aktive ettiği vs. baskıladığı yolaklar)
# yansıtabileceği gerçeğini göz ardı eder ve yön-spesifik sinyalleri
# sulandırabilir/maskeleyebilir. Bu yüzden GO ve KEGG zenginleştirmesini
# "up" ve "down" alt kümeleri için ayrı çalıştırıyoruz. Universe (arka
# plan) her iki yön için de aynı kalıyor: universe bir yön seçimi değil,
# "hangi genler test edilebilirdi" sorusunun cevabıdır.
sig_up_entrez <- clean_entrez_ids(
  sig_genes$gene_id[sig_genes$direction == "up"], "significant_genes.csv (up)"
)
sig_down_entrez <- clean_entrez_ids(
  sig_genes$gene_id[sig_genes$direction == "down"], "significant_genes.csv (down)"
)

cat(sprintf(
  "\nGeçerli Entrez ID sayısı: yukarı regüle = %d, aşağı regüle = %d, universe (tüm test edilenler) = %d\n",
  length(sig_up_entrez), length(sig_down_entrez), length(universe_entrez)
))

# "Sessizce geçme, uyar" kuralı: enrichGO/enrichKEGG'e anlamlı sayıda
# gen vermeden çağrı yapmak istatistiksel olarak neredeyse hiçbir zaman
# anlamlı bir zenginleştirme üretmez ve bu durumu sonuç tablosuna
# bakmadan fark etmek zor olur.
MIN_SIG_GENES <- 10
check_gene_count <- function(entrez_ids, label) {
  if (length(entrez_ids) < MIN_SIG_GENES) {
    cat(sprintf(
      "UYARI: %s için yalnızca %d anlamlı gen var (eşik: %d). Zenginleştirme sonuçları güvenilir olmayabilir.\n",
      label, length(entrez_ids), MIN_SIG_GENES
    ))
  }
  if (!all(entrez_ids %in% universe_entrez)) {
    n_outside <- sum(!(entrez_ids %in% universe_entrez))
    cat(sprintf(
      "UYARI: %s içindeki genlerden %d tanesi universe listesinde yok (beklenmiyordu).\n",
      label, n_outside
    ))
  }
}
check_gene_count(sig_up_entrez, "yukarı regüle genler")
check_gene_count(sig_down_entrez, "aşağı regüle genler")

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
# qvalueCutoff = 0.05: clusterProfiler'ın varsayılanı (0.2) yerine, üst
# akıştaki (02/03) anlamlılık eşiğiyle (padj < 0.05) tutarlı olacak
# şekilde bilinçli olarak 0.05 seçildi; aksi halde zenginleştirme adımı
# DGE adımından daha gevşek bir eşikle "anlamlı" terim üretip yanıltıcı
# olabilirdi.
run_go_enrichment <- function(gene_set, suffix, direction_label) {
  cat(sprintf("\n--- GO (BP) zenginleştirme çalışıyor: %s ---\n", direction_label))
  ego <- enrichGO(
    gene = gene_set,
    universe = universe_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE  # gen sembollerini de CSV'ye ekler, Entrez ID'den daha okunur
  )
  go_out <- file.path(results_dir, sprintf("go_enrichment_%s.csv", suffix))
  write.csv(as.data.frame(ego), file = go_out, row.names = FALSE)
  cat(sprintf("[tamam] GO sonuçları yazıldı: %s (%d terim)\n", go_out, nrow(as.data.frame(ego))))
  report_enrichment(ego, sprintf("GO terimi (%s)", direction_label))
  save_dotplot(
    ego,
    file.path(results_dir, sprintf("go_dotplot_%s.png", suffix)),
    sprintf("GO Biological Process zenginleştirme (%s)", direction_label)
  )
  ego
}

ego_up <- run_go_enrichment(sig_up_entrez, "up", "yukarı regüle")
ego_down <- run_go_enrichment(sig_down_entrez, "down", "aşağı regüle")

## 2) KEGG yolak zenginleştirme -----------------------------------------
# enrichKEGG, KEGG'in kendi REST API'sine canlı bir ağ isteği atar
# (organism="hsa" -> "Homo sapiens" KEGG kodu); internet yoksa veya
# KEGG servisi yanıt vermezse hata fırlatır. Bu hatayı sessizce yutmak
# yerine yakalayıp açıkça raporluyoruz, ama GO sonuçlarını (yukarıda
# zaten yazıldı) etkilememesi için script'in tamamını durdurmuyoruz.
run_kegg_enrichment <- function(gene_set, suffix, direction_label) {
  cat(sprintf(
    "\n--- KEGG zenginleştirme çalışıyor (%s, KEGG REST API'sine bağlanıyor) ---\n",
    direction_label
  ))
  ekegg <- tryCatch(
    enrichKEGG(
      gene = gene_set,
      universe = universe_entrez,
      organism = "hsa",
      keyType = "kegg",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05
    ),
    error = function(e) {
      cat(sprintf(
        "HATA: KEGG zenginleştirme (%s) başarısız oldu (muhtemelen ağ/KEGG API erişim sorunu): %s\n",
        direction_label, conditionMessage(e)
      ))
      NULL
    }
  )

  if (is.null(ekegg)) {
    cat(sprintf("UYARI: KEGG sonuçları (%s) üretilemedi.\n", direction_label))
    return(invisible(NULL))
  }
  ekegg_readable <- tryCatch(
    setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
    error = function(e) ekegg  # gen sembol dönüşümü başarısız olursa Entrez ID'lerle devam et
  )
  kegg_out <- file.path(results_dir, sprintf("kegg_enrichment_%s.csv", suffix))
  write.csv(as.data.frame(ekegg_readable), file = kegg_out, row.names = FALSE)
  cat(sprintf(
    "[tamam] KEGG sonuçları yazıldı: %s (%d yolak)\n",
    kegg_out, nrow(as.data.frame(ekegg_readable))
  ))
  report_enrichment(ekegg_readable, sprintf("KEGG yolağı (%s)", direction_label))
  save_dotplot(
    ekegg_readable,
    file.path(results_dir, sprintf("kegg_dotplot_%s.png", suffix)),
    sprintf("KEGG yolak zenginleştirme (%s)", direction_label)
  )
  ekegg_readable
}

ekegg_up <- run_kegg_enrichment(sig_up_entrez, "up", "yukarı regüle")
ekegg_down <- run_kegg_enrichment(sig_down_entrez, "down", "aşağı regüle")
