Normalization is the first important step required in analyzing RNA-seq in order to make valid comparisons across different samples.

This ultimately boils down to normalizing read count data. 
 - Do low read counts indicate a rare transcript? Or should they be discarded due to the uncertainty in their quantification.
 
 RNA-Seq and MicroArrays have the same problem  when it comes to detecting rare transcripts:
  - RNA-seq : Replicates of the same sample might not have a specific transcript if it is rare
  
There doesn't seem to be any standardized method for normalization so I have been testing as many different techniques as I can in a manner that makes at least a little bit of sense