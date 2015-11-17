(* config.ml *)



(*********************************************************)
(* Structure for the parameters of the overlap program   *)
(*********************************************************)
type 'a context_t =
	{ 
	  mutable file:string;			(* input bedpe file of the generalized exon to exon junctions with 24bp for donor and acceptor site *)
	  mutable outfile: string;              (* output bedpe file of the same exon to exon junctions with additional duplseq information *)
	} 



(********************************************************)
(* overlap context, these are the default parameters    *)
(********************************************************)
let context = 
  {	
    file = "";
    outfile = "";
  };;


let usage = 
"                   duplseq_gal - August 2014

Usage : "^(Filename.basename (Sys.argv.(0)))^" file [options] 

file is a bedpe file of generalized exon-exon junctions where the first block represents the 5' exonic segment 
until the donor site excluding it and where the second block represents the 3' exonic segment until the acceptor 
site excluding it.
Each row must contain at exactly 12 fields (1 denotes 5' exonic segment and 2 denotes 3' exonic segment):
- chr1
- start1
- end1
- chr2
- start2
- end2
- name
- score
- strand1
- strand2 
- the 24bp genomic sequence surrounding the donor site (12bp in exon, 12bp outside exon),
- the 24bp genomic sequence surrounding the acceptor site (12bp outside exon, 12bp in exon).

For each exon-exon junction of the input file, the program reports sequences of at least 5bp that are:
- present both in the donor and in the acceptor 24bp sequences, 
- either overlapping or adjacent to the 2bp junction of the 24bp sequence,
- not containing any bp present outside exons on both sides.

** file must be provided in bedpe format.
** [options] can be:
   -o outfile:   outfile is the name of the file where the output of duplseq_gal will be displayed.
                 -> default is file_withduplseq.bedpe
"

(***********************************************************************)
(* Read the arguments from the command line and updates the context    *)
(***********************************************************************)
let read_commandline () =
  let u = try 
      begin
	context.file <- Sys.argv.(1);
      end
    with
      | Invalid_argument s -> Common.print_error usage
  in
  
  (* we start reading the arguments from the 3rd one since the two first are compulsory and are the input files *)
  let argnum = ref 1 and ok = ref true in

  (* This function returns the next argument. mustbe says whether this argument must exist or not.
     In general mustbe is true since the arguments go by pairs *)
  let getarg mustbe =
    incr argnum; 
    if !argnum < (Array.length Sys.argv) then 
      Sys.argv.(!argnum)
    else if mustbe then 
      raise Not_found 
    else
      ""
  in
    (* Actually reading each of the arguments *)
    try 
      while (!ok) do
	match (getarg false) with
	  | "-o" -> context.outfile <- getarg true
	  | "-h"
	  | "--help" -> Printf.printf "%s\n" usage; exit 1; 
	  | ""	-> ok := false
	  | s	-> failwith s
      done;
      Common.print_log "Command line is read\n";
    with
      | Not_found -> Common.print_error ("Missing parameter at the end of the command line\n"^usage^"\n");
      | Failure "int_of_string" -> Common.print_error ("Syntax error (incorrect integer value)\n"^usage^"\n");
      | Failure s -> Common.print_error ("Syntax error ("^s^")\n"^usage^"\n");;


