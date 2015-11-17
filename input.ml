(* input.ml *)
(* reads input file and store data in appropriate structures *)


open Common
open Feature



(* record_of_line_bedpe takes as input a line of a bedpe file as a list of strings split by the "\t" separator, 
   and outputs a feature object. 
   Note that the bedpe file must have exactly 12 fields where the last 2 correspond to two 24bp genomic sequence
   at the donor and acceptor sites of the generalized exon to exon junctions
*)
let record_of_line_bedpe lline =
  let tab = Array.of_list lline in
  let n = Array.length tab in
    if (n != 12) then
      begin
	failwith "record_of_line_bedpe: Syntax error in your bedpe file: it should contain exactly 12 fields separated by tabs"
      end
    else
      Feature.create 
	tab.(0) 
	(int_of_string tab.(1))
	(int_of_string tab.(2))
	tab.(3) 
	(int_of_string tab.(4)) 
	(int_of_string tab.(5)) 
	tab.(6)
	tab.(7)
	(strand_of_string tab.(8)) 
	(strand_of_string tab.(9)) 
	[tab.(10) ; tab.(11)] 
	

(* This function reads the lines of a bedpe file and stores them in a list of features. 
   infile is the name of the input bedpe file. *)
let make_feat_list infile = 
  let inchan = open_in infile and lfeat= ref [] and stop = ref false in
    while (not !stop) do
	(
	  try
	    let currline = (split '\t' (input_line inchan)) in
	    lfeat:=(record_of_line_bedpe currline)::(!lfeat) 
	  with
	    | End_of_file -> stop:=true
	)
      done;
      !lfeat









