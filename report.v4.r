
argv <- commandArgs(T)
if (length(argv)!=2) {
		stop("\n\t\tUsage: Rscript init.r [profile] [phenotype]\n\n",call. =F)
}
library('shiny')
library('DT')  # for table section
library('vegan')
library('randomForest')
library('pROC')
library('heatmap3')
library('ade4')
#library('fpc')
#library('e1071')  #SVM
#library('rpart')  #CART
#library('adabag') #adaBoosting
#library('ggplot2')
#library('glmnet') #feature selection by lasso
#library('RSNNS')  #neutral network
#library('hclust')
#library('pam')
#library('plyr')
library('psych')
library('igraph')
library('fdrtool')

prof <- as.matrix(read.table(argv[1],header=T,row.names=1));  t.prof <- t(prof)
phenotype <- read.table(argv[2],header=T,row.names=1)
p <- phenotype[,1]
uni.p <- unique(p)
a <- which(p==uni.p[1])
b <- which(p==uni.p[2])

num.p <- sapply(1:ncol(phenotype),function(x) as.numeric(phenotype[,x]))
colnames(num.p) <- colnames(phenotype)
rda.scaling.options <- c("关注物种","关注样品","同时关注物种与样本")

corMethod <- c("spearman","pearson","kendall")


# ~~~~~~~ wilcox rank sum test ~~~~~~
rst <- function(prof,a,b,uni.p){
	pval <- sapply(1:nrow(prof),function(x) wilcox.test(prof[x,a],prof[x,b],alternative = "two.sided")$p.value)
	pval[is.na(pval)] <- 1
	fdr <- fdrtool(pval,statistic="pvalue")$lfdr
	mean.A <- rowMeans(prof[,a])
	mean.B <- rowMeans(prof[,b])
	freq.A <- sapply(1:nrow(prof),function(x) length(which(prof[x,a]>0))/ncol(prof[,a]))
	freq.B <- sapply(1:nrow(prof),function(x) length(which(prof[x,b]>0))/ncol(prof[,b]))
	enrichment <- as.character(sapply(1:nrow(prof),function(x) if(mean.A[x] > mean.B[x]) uni.p[1] else uni.p[2]))
	mat <- cbind(pval,fdr,mean.A,mean.B,freq.A,freq.B)
	rownames(mat) <- rownames(prof)
	mat <- mat[order(pval),]
	mat <- as.data.frame(mat)
	mat <- cbind(enrichment,mat)
	return(mat)
}

cor.sig <- function(prof,method,cor.r,cor.p){
	cor.obj <- corr.test(prof,method=method)
	cor.mat <- cor.obj$r
	pval.mat <- cor.obj$p
	
	pval.mat[pval.mat <= cor.p] <- -1
	pval.mat[pval.mat > cor.p] <- 0
	pval.mat <- abs(pval.mat)

	cor.mat[abs(cor.mat) >= cor.r] <- -1
	cor.mat[abs(cor.mat) < cor.r] <- 0
	cor.mat <- abs(cor.mat)

	union.mat <- pval.mat * cor.mat
	diag(union.mat) <- 0
	return(union.mat)
}
#~~~~~~~~~~~ classificatoin~~~~~~~~~~~~~
# ~~~random forest
rf_model <- function(prof,p,ntree,nPerm){
	rf.prof <- t(as.data.frame(prof))
	colnames(rf.prof) <- gsub('-|;|[[]|[]]','.',rownames(prof))
	rf.prof <- cbind(p,rf.prof)
	colnames(rf.prof)[1] <- "group"
	set.seed(211)
	fit <- randomForest(factor(group)~., data = rf.prof, 
				ntree=ntree,proximity=T,importance=T,nPerm=nPerm)
	return(fit)
}


cls.model <- c("Random Forest","SVM","naive Bayesian","adaBoosting","CART","Neutral Network")

heatmap.col <- list(colorRampPalette(c("navy", "white","firebrick3"))(1024),
				colorRampPalette(c("white","black"),space="Lab")(256),
				colorRampPalette(c("black","green","green"),space="Lab")(125),
				colorRampPalette(c("#009933","white","#CC9933"),space="Lab")(256),
				topo.colors(256),heat.colors(256),cm.colors(256)
				)
names(heatmap.col) <- c('default','white-black','black-green','green-yellow','topo','heat','cm')

# ~~~ end classification ~~~~~~~~~~~~

pca <- function(prof){
	res <- dudi.pca(prof,scan=F,center=T,scale=F,nf=2)
	return(res)
}


# build a shiny app object which can be used by runApp()
app <- shinyApp(
	server = shinyServer(function(input, output) {

  		output$distPlot <- renderPlot({
			i <- input$densPlotIndex
		    x <- as.numeric(prof[i,])
			bins <- seq(min(x), max(x), length.out = input$bins + 1)
	    	hist(x, breaks = bins, col = 'skyblue', border = 'white',
				main = rownames(prof)[i]
			)
  		})

		output$personPie1 <- renderPlot({
			cn <- input$personIndex1
			value <- prof[,cn]
			names(value) <- gsub('.*g__','',names(value))
			pie(sort(value,decreasing=T)[1:10])
		})
		output$personPie2 <- renderPlot({
			cn <- input$personIndex2
			value <- prof[,cn]
			names(value) <- gsub('.*g__','',names(value))
			pie(sort(value,decreasing=T)[1:10])
		})
		output$personPie3 <- renderPlot({
			cn <- input$personIndex3
			value <- prof[,cn]
			names(value) <- gsub('.*g__','',names(value))
			pie(sort(value,decreasing=T)[1:10])
		#	title(colnames(prof)[cn],line=-1)
		})
  		
  		
		output$boxPlot <- renderPlot({
		    A <- vegan::diversity(t(prof[,a]))
		    B <- vegan::diversity(t(prof[,b]))
			
			par(mar=c(5,5,4,2))
	    	boxplot(list(A,B), col = c('skyblue','orange'),
				names=uni.p, ylab = "shanon diversity")
			
  		})
  		
		output$densPlot <- renderPlot({
			i <- input$densPlotIndex
		    x <- as.numeric(prof[i,])

			den1 <- density(x[a])
			den2 <- density(x[b],bw=den1$bw)
			par(mar=c(5,2,4,1))
	    	plot(den1,main = rownames(prof)[i],ylim=c(0,max(den1$y)*1.2) )
			lines(den2,col='skyblue')

  		})

## ~~~~~  reactive object ~~~~~~~~		
		rf <- reactive({
			rf_model(prof,p,input$rf.ntree,input$rf.nPerm)
		})
		PCA <- reactive({
			pca(t.prof) 
		})



## ~~~~~  reactive object ends here
		output$randomForestErRate <- renderPlot({
			plot(rf())
		})
		output$randomForestImp <- renderPlot({
			imp <- importance(rf())
			i.gini <- 4
			imp <- as.data.frame(imp[order(-imp[,i.gini]),])
	        imp2 <- imp[1:20,i.gini]
			name <- rownames(imp)[1:20]			
			name <- gsub('.*f__','',name)

			par(mar=c(5,18,4,2))	
			barplot(as.numeric(imp2),col='grey80',horiz=T,las=2,
				xlab = colnames(imp)[i.gini],
				names = name 	
				)
		})
		output$randomForestROC1 <- renderPlot({
			fit <- rf_model(prof,p,input$rf.ntree,input$rf.nPerm)
		    pred<-predict(fit,type="prob")
		    p_true<-as.numeric(factor(p))
		    p_pred<-as.numeric(pred[,2])
		    roc(p_true,p_pred,plot=T,print.thres=T,print.auc=T)
		})
		
		output$randomForestROC2 <- renderPlot({
			fit <- rf_model(prof,p,input$rf.ntree,input$rf.nPerm)
		    pred<-predict(fit,type="prob")
		    p_true<-as.numeric(factor(p))
		    p_pred<-as.numeric(pred[,2])
			roc1 <- roc(p_true,p_pred,percent=T,plot=T,print.auc=T)
			sens.ci <- ci.se(roc1)
			plot(sens.ci,type="shape",col="lightgrey",main=ci.auc(roc1))
		})

		output$pcaPlot <- renderPlot({
			pca.res <- PCA() 
			pca.li <- pca.res$li[,1:2]
			pca.eig <- pca.res$eig[1:20]

			pch <- rep(17,length(p))
			pch[b] <- rep(15)
			par(mfrow=c(1,2))
			plot(pca.li,col=pch-13,pch=pch,xlab="PC1",ylab="PC2",main="")
			plot(pca.eig,type='b',col='red')
		})

		output$svdPlot <- renderPlot({
			svd.res <- svd(t.prof)
			svd.u <- svd.res$u[,1:2]
			svd.d <- svd.res$d[1:20]
			pch <- rep(17,length(p))
			pch[b] <- rep(15)
			par(mfrow=c(1,2))
			plot(svd.u,col=pch-13,pch=pch,xlab="Axis1",ylab="Axis2",main="")
			plot(svd.d,col='orange',type='b')
		})

		output$rdaPlot <- renderPlot({
			rd <- rda(t.prof,num.p)
			sc.i <- which(input$rda.scaling==rda.scaling.options)
			plot(rd,scaling=sc.i,type='n',cex.axis=1,cex.lab=1,main="RDA diagram")
			points(rd,scaling=sc.i,display="si",col="royalblue",pch=20,axis.bp=F,cex=1)
			points(rd,scaling=sc.i,display="sp",col="red",pch=17,axis.bp=F,cex=1)
			text(rd,display="bp",col='black',cex=1,axis.bp=F)
			if(input$showText.sp==T) 
				text(rd,display="sp",col='red',cex=0.5,axis.bp=F)
			if(input$showText.si==T) 
				text(rd,display="si",pos=1,offset=0.5,col='royalblue',cex=0.5,axis.bp=F)
			
		})		

		output$rdaUI <- renderUI({
			plotOutput('rdaPlot',width=paste(input$rda.width,'px',sep=''),height=paste(input$rda.height,'px',sep=''))
		})
		
		rank.test <- reactive({
			rst(prof,a,b,uni.p)
		})	

		output$dynamic.heatmap <- renderUI({
			plotOutput('heatmap',width=paste(input$heatmap.width,'px',sep=''),height=paste(input$heatmap.height,'px',sep=''))
		})

		output$heatmap <- renderPlot({
			rt.mat <- rank.test()
			x <- prof[rt.mat[,2]<0.01 & rt.mat[,3]<0.05,]
			rownames(x) <- gsub('.*f__','',rownames(x))

			c.index <- which(names(heatmap.col)==input$heat.col)
			heatmap3(x,Colv=NA,
				col = heatmap.col[[c.index]],
				cexRow=1,cexCol=1,
				ColSideLabs=NA,
				labCol="",scale='row',
				ColSideColors=as.numeric(factor(p)),
				margin=c(1,8)				
			)

		})

		output$correlationPlot <- renderPlot({
			marker.prof <- t.prof[,sample(1:ncol(t.prof),50)]
			colnames(marker.prof) <- paste(gsub('.*f__','',colnames(marker.prof)),1:ncol(marker.prof),sep='.')
			res.mat <- cor.sig(marker.prof,input$cor.method,input$cor.r,input$cor.p)
			graph <- graph.adjacency(res.mat,mode='undirected')

			fc <- cluster_fast_greedy(graph)
			group <- communities(fc)
			plot(graph,
				layout = layout_with_fr,
				vertex.size=2,
				vertex.color = 'green',
				vertex.label.dist = 0.15,
				vertex.label.cex = 1,
				vertex.label.degree = 90,
				edge.width = 0.7,
				mark.groups = group[sapply(1:length(group),function(x) length(group[[x]])>2)]
			)
		
		})

		output$sampleBarplot <- renderPlot({
			top.n <- input$bar.topN
	
		    index <- order(-rowMeans(prof))
			i.top <- index[1:top.n]
			i.other <- index[(top.n+1):length(index)] 
		    x <- prof[i.top,]
			x <- rbind(x,colSums(prof[i.other,]))
			set.seed(88520)
			color <- sample(rainbow(top.n+1),top.n+1)
			
			ai <- sample(a,input$bar.aN)
			bi <- sample(b,input$bar.bN)
			x2 <- cbind(x[,ai],x[,bi])
			barplot(x2,col = color,las=2)	
#			legend("topleft",legend=rownames(prof)[index[1:top.n]],pch=15,col=color)
		})

# ~~~~  report download session starts here  ~~~~~~~~~~

		regFormula <- reactive({
		    as.formula('mpg ~ cyl')
		})


	  	output$downloadReport <- downloadHandler(
    		filename = function() {
      			paste('Microbiota_analysis_report',input$version,sep='.', switch(
        		input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
      			))
    		},

    	content = function(file) {
      		src <- normalizePath('report.Rmd')

	      # temporarily switch to the temp dir, in case you do not have write
    	  # permission to the current working directory
      		owd <- setwd(tempdir())
	      	on.exit(setwd(owd))
    	  	file.copy(src, 'report.Rmd', overwrite = TRUE)

      		library(rmarkdown)
      		out <- rmarkdown::render('report.Rmd', 
						switch(input$format,
				        	PDF = pdf_document(), HTML = html_document(), Word = word_document())
#      					params = list(PCA=pcaPlot(),SVD=svdPlot())	
			)

      		file.rename(out, file)
    	}
	)	


# ~~~ report download session ends here

# ~~~  table download and upload ~~~~
  datasetInput <- reactive({
    switch(input$download.tb,
           "phenotype" = phenotype,
           "profile" = prof,
           "classification_importance" = importance(rf()),
           "rank_sum_test" = rank.test())
})
		output$downloadData <- downloadHandler(
    		filename = function() { 
		 	paste(input$download.tb, 'csv', sep='.') 
	 		},
    		content = function(file) {
		      write.csv(datasetInput(), file)
    		}
		)	


# ~~~~~~~  clustering by kmeans with two demension variables
	 selectedData <- reactive({
    	t.prof[,c(input$xcol, input$ycol)]
  		})

  	clusters <- reactive({
    	kmeans(selectedData(), input$clusters)
  	})

	  output$clusterPlot <- renderPlot({
    	palette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
     	 "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"))

    	par(mar = c(5.1, 4.1, 0, 1))
	    plot(selectedData(),
         	col = clusters()$cluster,
        	 pch = 20, cex = 3)
    	points(clusters()$centers, pch = 4, cex = 4, lwd = 4)
  	})	


# ~~~~~~~~~~ table section start ~~~~~~~~~~~~~~
	rank.test <- reactive({
		rst(prof,a,b,uni.p)
	})	

  output$phenoTbl <- DT::renderDataTable({
    inFile <- input$mappingTbl
    if (is.null(inFile))
		return(NULL)
    DT::datatable(read.csv(inFile$datapath, header=input$headerMapping, sep=input$sepMapping,quote=input$quoteMapping),
    	options = list(orderClasses = TRUE))
  })

  output$profTbl <- DT::renderDataTable({
    inFile <- input$profileTbl
    if (is.null(inFile))
		return(NULL)
    DT::datatable(read.csv(inFile$datapath, header=input$headerProfile, sep=input$sepProfile,quote=input$quoteProfile),
    	options = list(orderClasses = TRUE))
  })

  # customize the length drop-down menu; display 5 rows per page by default
  output$rtTbl <- DT::renderDataTable({
    DT::datatable(rank.test(), options = list(lengthMenu = c(10, 30, 50), pageLength = 10))
  })

  output$impTable <- DT::renderDataTable({
	imp <- importance(rf())
	imp <- imp[order(-imp[,4]),3:4]
    DT::datatable(imp, options = list(lengthMenu = c(10, 30, 50), pageLength = 5))
  })


# ~~~~~~~~ table section end ~~~~~~~~


}),  # server function end here

	# Define UI for application that draws a histogram
	ui = shinyUI(fluidPage(
	titlePanel("谱元生信分析平台v1.0"),
	
	  # Sidebar with a slider input for the number of bins
 	sidebarLayout(position="left",
    sidebarPanel(
			
		h3("堆砌柱形图"),
		fluidRow(
			column(4,sliderInput("bar.topN","丰度前几物种",
                	  min = 1,max = 20,value = 8)),
			column(4,sliderInput("bar.aN","A组显示样本数",
                	  min = 1,max = length(a),value = 10)),
			column(4,sliderInput("bar.bN","B组显示样本数",
                	  min = 1,max = length(b),value = 10))
		),

		h3("热图参数"),
		fluidRow(
			column(5,sliderInput("bins","数值转换",
                	  min = 0 ,max = 50,value = 20)),
			column(5,sliderInput("b2","数值转换",
                	  min = 0 ,max = 50,value = 20))
		),
			selectInput('heat.col',"热图颜色方案",names(heatmap.col)),
			numericInput('heatmap.width',"热图宽度",'800'),
			numericInput('heatmap.height',"热图长度",'600'),

		h3("特定物种分组比较"),
    	selectInput('densPlotIndex', '物种选择', rownames(prof),selected=rownames(prof)[2]),
		fluidRow(
    		column(4,selectInput('personIndex1', '个体A选择', colnames(prof),selected=colnames(prof)[1])),
    		column(4,selectInput('personIndex2', '个体B选择', colnames(prof),selected=colnames(prof)[2])),
    		column(4,selectInput('personIndex3', '个体C选择', colnames(prof),selected=colnames(prof)[3]))
		),
   		h3('分类器参数'),
		
	    selectInput('cls.model', '分类模型', cls.model,selected=cls.model[1]),
		numericInput('rf.ntree','随机森林决策树个数','500'),
		numericInput('rf.nPerm','随机森林迭代次数(nPerm)','1'),

		h3("RDA与CCA分析"),
		fluidRow(
			column(4,selectInput('rda.scaling',"scaling",rda.scaling.options,selected=rda.scaling.options[[1]])),
			column(4,checkboxInput('showText.sp', '显示物种名', FALSE)),
			column(4,checkboxInput('showText.si', '显示样本名', FALSE))
		),
		fluidRow(
			column(5,numericInput('rda.height',"RDA图高度","600")),
			column(5,numericInput('rda.width', 'RDA图宽度', "600"))
		),

   		h3('数据上传与下载'),
		fluidRow(
			column(5, selectInput("download.tb","选择要下载的数据表",choices = c("classification_importance","rank_sum_test","profile","phenotype"))),
			br(),
			column(4, downloadButton("downloadData","下载"))
		),
		fluidRow(
			column(5, 
				fileInput("profileTbl", label = h3("丰度表上传"),accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
				checkboxInput('headerProfile', 'Header', TRUE),
				radioButtons('sepProfile', 'Separator',c(Comma=',',Semicolon=';',Tab='\t'),','),
				radioButtons('quoteProfile', 'Quote',c(None='','Double Quote'='"','Single Quote'="'"),'"')
			),
			column(5, 
				fileInput("mappingTbl", label = h3("表型文件上传"),accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
				checkboxInput('headerMapping', 'Header', TRUE),
				radioButtons('sepMapping', 'Separator',c(Comma=',',Semicolon=';',Tab='\t'),','),
				radioButtons('quoteMapping', 'Quote',c(None='','Double Quote'='"','Single Quote'="'"),'"')
			)
		),

  		h3('报告下载'),
      		helpText("请将报告下载为word文档格式,PDF格式开发中"),
      		radioButtons('format', 'Document format', c('Word', 'HTML', 'PDF'),inline = TRUE),
			textInput('version',"报告版本号","v1"),
			downloadButton('downloadReport'),

		br(),
		br(),
		br(),
		br(),
		br(),
   		h3('聚类参数'),
    	selectInput('xcol', 'X Variable', rownames(prof)),
	    selectInput('ycol', 'Y Variable', rownames(prof),selected=rownames(prof)[[2]]),

	    fluidRow(
			column(3,numericInput('clusters', label = h4('聚类簇数'), 3,min = 1, max = 9)),

			column(3,radioButtons("clusterModel",label = h4("聚类模型"),
	 	         choices = list("K-means" = 1,"PAM" = 2),
 		         selected = 1)),
			column(3,radioButtons("clusterDist",label = h4("距离度量"),
	 	         choices = list("Euclidean" = 1,"Bray-curtis" = 2),
 		         selected = 1)),
			column(3,br(),
				actionButton("action", label = "重置"),
				br(),br(),  #按键上下间距
			submitButton("确定"))
		),
		
		h3("特征菌群相关性网络"),
		fluidRow(
			column(4,selectInput("cor.method","相关性方法",corMethod)),
			column(4,numericInput("cor.p","p value",0.001)),
			column(4,numericInput("cor.r","相关系数",0.4))
		)
	), # sideBarPanel ends here

    # Show a plot of the generated distribution
    mainPanel(
		h2("1.整体菌群分析"),
		h3("1.1 群落堆砌柱形图"),
			plotOutput("sampleBarplot"),
		h5("barPlot legend:"),
		h5(rownames(prof)[1]),
		h5(rownames(prof)[2]),
		fluidRow(
			column(4,plotOutput("personPie1")),
			column(4,plotOutput("personPie2")),
			column(4,plotOutput("personPie3"))
		),
		h3("1.2 PCA主成分分析"),
			plotOutput("pcaPlot",width="100%"),	
			includeText("txt/pca.txt"),	
			plotOutput("svdPlot",width="100%"),	
			includeText("txt/svd.txt"),	
		h2("2.菌群分组比较"),
		h3("2.1 香农多样性指数"),
		fluidRow(
			column(4,plotOutput("boxPlot")),
			column(8,plotOutput("densPlot"))
		),
			plotOutput("distPlot"),
# ~~~~~~~~~~~~ This section arrange tables for raw data,rank test
		h3("2.2 分组差异检验"),
	fluidPage(
	#  title = 'DataTables title',
    	mainPanel(
      	tabsetPanel(
        	id = 'dataset',
        	tabPanel('物种分类预测权重', DT::dataTableOutput('impTable')),
        	tabPanel('phenotype', DT::dataTableOutput('phenoTbl')),
        	tabPanel('profile', DT::dataTableOutput('profTbl')),
        	tabPanel('wilcox test', DT::dataTableOutput('rtTbl'))
      	)
    	)
	),

# ~~~~~~~~~~~ end table section ~~~~~~~~~~~~~~~~~
		h3("差异物种热图可视化"),
			uiOutput("dynamic.heatmap"),
		h3("RDA与CCA分析"),
			uiOutput("rdaUI"),
			h4("上图中红点表示物种，蓝点表示样本，箭头方向与样本、物种的相对位置表示各因素间相关性强弱"),
			includeText("txt/cca.txt"),

		h2("3. 疾病分类器构建"),
			plotOutput("randomForestErRate"),
			fluidRow(				
				column(5,plotOutput("randomForestROC1")),				
				column(5,plotOutput("randomForestROC2"))
			),				
			plotOutput("randomForestImp"),				

		h2("4. 探索性分析"),
		h3("4.1 无监督聚类"),
			plotOutput("clusterPlot"),
		h3("4.2 菌群特征选择"),
		h3("5. 特征菌群相关性网络"),
			plotOutput("correlationPlot",height="800px",width="800px"),
		br()
    )
  	)
	))
)


runApp(app,host='192.168.1.60',launch.browser=F)



