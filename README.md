# R.Codes
#Heatmap question
چطور فقط یک لیبل برای سمپل ها داشته باشم یعنی فقط ERStatus در نهایت در هیت مپ نشون بده و بدون اینکه خوشه بندی بشه یعنی تمام ERpos ها سمت راست و همه ERneg ها سمت چپ؟ بدون اینکه mislabeling اتفاق بیفته من این بخش قبل رسم هیت مپ اضافه کردم ولی همچنان هر دو لیبل ER+ و ER- نشون میده وقتی بخوام حذفش کنم خطا میده.
*Step 1: Identify samples with ER+ = 1
ER_positive_samples <- rownames(metabric_annot_disc)[metabric_annot_disc$`ER+` == 1]
ER_negative_samples <- rownames(metabric_annot_disc)[metabric_annot_disc$`ER+` == 0]

*Step 2: Order columns in metabric_regact_disc with ER+ samples on the left
ordered_sample_names <- c(ER_positive_samples, ER_negative_samples)
metabric_regact_disc_ordered <- metabric_regact_disc$differential[ordered_sample_names,]

ER_annotation <- metabric_annot_disc[,c("ER+","ER-")]
ER_annotation_colors <- list(
  "ER+" = c("0" = "lightgrey", "1" = "blue"),
  "ER-" = c("0" = "lightgrey", "1" = "red")
  )
# Experience
Always Create the code with a distinct name; instead of changing the code for every situation; just create a new project in the name that points out the corresponding label
