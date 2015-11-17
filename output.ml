(* output.ml *)


open Feature


(* print takes as input:
   - o which is the output channel corresponding to the file where we need to write file1 with additional info 
   - flist which is the feature list corresponding to the features of the input file with added information 
   and prints in o the features (however there is no order in the sequences printed out)
*)
let print o flist =
  List.iter (Feature.print_bedpe o) flist
  


 
