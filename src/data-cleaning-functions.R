#
# Title: Data cleaning functions
# Created: March 28th, 2022
# Last Updated: March 28th, 2022
# Author: Brandon Allen
# Objectives: Functions required to summarize the landscape long-form data
# Keywords: Landscape summaries
#

#######################
# Landscape summaries # 
#######################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

landscape_hf_summary <- function(data.in, landscape.lookup, class.in, class.out) {
        
        # Matching of lookup tables and merging native features
        landscape.lookup <- landscape.lookup[landscape.lookup[, class.in] %in% colnames(data.in), ]
        landscape.clean <- matrix(nrow = nrow(data.in), ncol = length(unique(landscape.lookup[, class.out])))
        
        for (abmi.coef in 1:length(unique(landscape.lookup[, class.out]))) {
                
                coef.temp <- as.character(landscape.lookup[landscape.lookup[, class.out] %in% as.character(unique(landscape.lookup[, class.out]))[abmi.coef], class.in])
                
                if(length(coef.temp) == 1) {
                        
                        landscape.clean[, abmi.coef] <- data.in[, coef.temp]
                        
                } else {
                        
                        landscape.clean[, abmi.coef] <- rowSums(data.in[, coef.temp])
                        
                }
                
        }
        
        colnames(landscape.clean) <- as.character(unique(landscape.lookup[, class.out]))
        rownames(landscape.clean) <- rownames(data.in)
        
        return(landscape.clean)
        
}
