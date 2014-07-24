test1 <- c()
count = 1
para <- for(i in affil){
                for(j in cou){
                        if(grepl(j,i) == TRUE){
                                count <- count + 1
                                test1 <- paste(test1,j)
                        }
                }
}

