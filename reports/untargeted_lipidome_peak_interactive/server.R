library(shiny)



function(input, output, session) {
        output$hist <- renderPlot({
                hist(lipid$fdata$qc_cv, main = "The distribution of CVs")
        })
        output$cv30 <- renderText({
                sum(lipid$fdata$qc_cv>30)
        })
        output$blank <- renderText({
                sum(lipid$fdata$sample_mean < lipid$fdata$blank_min)
        })
        # output$lipids <- renderUI({
        #         selectInput("feature", "Lipid species", choices = lipidNames)
        # })
        output$pie <- renderPlot({
                plotPie()
        })
        output$bars <- renderPlot({
                feature <- featureNames(lipid2)[lipid2$fdata$name==input$feature]
                plotBars(feature)
        })
}
