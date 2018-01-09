legend("topright", legend=c("Tumor", "Normal"), fill=c("red", "blue"), bty="n")
par(mar=c(10,6,3,3))

barplot(as.numeric(c(b9.data["A_25_P00010504",],b9.data.normal["A_25_P00010504",])),las=3, ylab="Log2(Expression Level)", main="Expression of hsa-let-7f-5p\n(fold change=58)",col=c(rep("red",39), rep("blue",8)), names.arg=c(colnames(b9.data), colnames(b9.data.normal)))
