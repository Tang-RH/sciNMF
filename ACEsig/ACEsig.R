# This function is used to calculate the ACE score and classify the patients into high and low risk groups based on the median of the ACE score
# The input is expression matrix, whose rows represent genes and columns represent patients
# return value is a dataframe object, frist column is 'ACEScore' cacluated by the ACE model, the second column is 'Risk' classification 
ACEsig = function(mat){
    model = c('CALCRL'=0.516705634819962,
              'CST3'=-0.135377434640815,
              'DDIT4'=0.147011776933904,
              'GABARAP'=-0.503454274985267,
              'HOPX'=0.110456507784065,
              'ID2'=0.139516214077651,
              'IL32'=-0.310282298890454,
              'ITM2A'=-0.151790695494774,
              'LSP1'=0.223798555072212,
              'NPDC1'=0.110819689132806,
              'SPINT2'=0.245331132788398,
              'SRGN'=0.382375238214482)
    genes = names(model)
    idx_g = genes %in% rownames(mat)
    if(any(!idx_g)){
        warning('Gene: ', paste0(genes[!idx_g],collapse = ', '), ' are(is) absent in your expression matrix')
    }
    coeff = model[idx_g]
    Score = apply(mat[names(coeff),], 2, function(x){
        sum(x*coeff)
    })
    cat('Median Cutoff is ', median(Score))
    Risk = ifelse(Score > median(Score), 'High','Low')
    Patient = colnames(mat)
    df_res = data.frame(Patient = Patient, Score = Score, Risk = Risk)
    rownames(df_res) = Patient
    return(df_res)
}