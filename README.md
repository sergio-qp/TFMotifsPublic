# TFMotifs

## A bioinformatics pipeline for the delineation of co-binding transcription factor conformations

Every specific trait and function the body has is managed by the carefully regulated expression of 
genes. However, understanding how exactly this regulation is achieved is difficult since it is carried 
out by highly interconnected and complex systems of proteins called “transcription factors” (TFs). 
These proteins function by binding to specific DNA sequences called “motifs” on “non-coding” 
regions of the genome, sections that do not directly encode proteins, but rather act as regulatory 
spaces, determining to what extent any given gene is expressed. What must be noted here is that 
the binding of a TF to a non-coding region rarely directly leads to a simple activation or 
inactivation, rather several different TFs interact with each other in specific conformations, fine 
tuning the extent to which a gene is expressed at any given stage. The problem lies in that many 
complex diseases, such as type 2 diabetes (T2D), are caused by variations in the non-coding DNA 
sequences these TFs bind to. This makes them difficult to accurately describe and understand, as 
the interconnected nature and high sensitivity of TF systems makes them very time consuming to 
thoroughly examine. In our study, we propose an analysis pipeline through which we can easily 
highlight areas of interest within the range of possible interactions between cooperating TFs. 
Given that any one TF tends to interact with any other in a specific manner, we designed a 
computational method to systematically describe the physical arrangements of any given TF pair 
that happen more frequently than would be expected, suggesting some sort of preference. This 
method combines JASPAR motif data for our selected transcription factors with ChIP-seq peaks from
the Pasquali et al. paper "Pancreatic islet enhancer clusters enriched in type 2 diabetes
risk-associated variants". The pipeline is illustrated below.

<img width="757" height="840" alt="image" src="https://github.com/user-attachments/assets/dcfe758b-874f-4a77-922b-49873e361458" />

