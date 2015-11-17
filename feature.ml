(* feature.ml *)

type strand = | Forward  
              | Reverse
	      | Unstranded
	      | Unknown


let strand_of_string s =
  match s with
    | "++"
    | "+" -> Forward
    | "+-"
    | "-" -> Reverse
    | "." -> Unstranded
    | "?" -> Unknown
    |  _  -> failwith "strand_of_string : the strand should be provided as the 7th field of your gtf file\n" 


let strand_to_string = function
  | Forward -> "+"
  | Reverse -> "-"
  | Unstranded -> "."
  | Unknown -> "?"


let strand_to_string2 = function
  | Forward -> "p"
  | Reverse -> "m"
  | Unstranded -> "."
  | Unknown -> "?"




(* Corresponding to a bedpe format with 12 fields. 
   1 denotes the donor block while 2 denotes the acceptor block.
*)
module Feature =
struct
  type t = {
    seq1: string;                     
    start1: int;                        
    end1: int;           (* coord of 1st block = donor, 0-based *)               
    seq2: string;                     
    start2: int;                   
    end2: int;           (* coord of 2nd block = acceptor, 0-based *)        
    name: string;
    score: string;                  
    strand1: strand;
    strand2: strand;
    attlist: string list  (* the attribute list but where the two first attributes correspond to the donor and acceptor 24bp genomic sequences *)
  }

  let create se1 st1 en1 se2 st2 en2 nm sc str1 str2 al = 
    {
      seq1=se1; 
      start1=st1; 
      end1=en1; 
      seq2=se2; 
      start2=st2; 
      end2=en2; 
      name=nm; 
      score=sc; 
      strand1=str1; 
      strand2=str2; 
      attlist=al
    }
  
  let null = create "." (-1) (-1) "." (-1) (-1) "null" "0" Unstranded Unstranded []
    
  let setattlist al f = {f with attlist=al}

  (* from a feature attribute list makes a string for bedpe output *)
  let string_attlist_in_bedpe f = 
    let s = ref "" in
    let u = List.iter (fun v -> s:=(!s)^(" ")^(if(v="") then "." else v)) f.attlist in
      (!s)

   (* 
      print_bedpe takes as input:
      - o which is the output channel corresponding to the file where we need to write file with duplseq info
      - f which is a feature
      and prints in o the feature in bedpe format.
   *)
  let print_bedpe o f =
    Printf.fprintf o "%s\t%i\t%i\t%s\t%i\t%i\t%s\t%s\t%s\t%s\t%s\n" f.seq1 f.start1 f.end1 f.seq2 f.start2 f.end2 f.name f.score (strand_to_string f.strand1) (strand_to_string f.strand2) (string_attlist_in_bedpe f)
 
  let seq1 f = f.seq1
  let start1 f = f.start1
  let end1 f = f.end1
  let seq2 f = f.seq2
  let start2 f = f.start2
  let end2 f = f.end2
  let name f = f.name
  let score f = f.score
  let strand1 f = f.strand1
  let strand2 f = f.strand2
  let attlist f = f.attlist
end





module Exon =
struct
  type t = {
	chrom: string;
	gbeg: int;
	gend: int;
	str: strand;
	score: float;
	cat: string;
	trid: string;
	gnid: string;
  }

  let create chr gb ge st sc ct tr gn = { chrom=chr; gbeg=gb; gend=ge; str=st; score= sc; cat=ct; trid=tr; gnid=gn}
  let null = create "" (-1) (-1) Forward  0.0 "" "" ""

  (* when we compare two exons we compare them by gene first *)
  let compare e1 e2 = 
    if ((Pervasives.compare e1.gnid e2.gnid) !=0) then
      Pervasives.compare e1.gnid e2.gnid
    else
      begin
	if ((Pervasives.compare e1.gbeg e2.gbeg) !=0) then
	  Pervasives.compare e1.gbeg e2.gbeg
	else
	  Pervasives.compare e1.gend e2.gend
      end
    
  let chrom e = e.chrom
  let gbeg e = e.gbeg
  let gend e = e.gend
  let str e = e.str
  let score e = e.score
  let cat e = e.cat
  let trid e = e.trid
  let gnid e = e.gnid
end 


module Transcript =
struct
  type t = {
	chrom: string;
	gbeg: int;
	gend: int;
	str: strand;
	cat: string;
	exlist: Exon.t list;
	trid: string;
	gnid: string;
  }
	  
  let create chr gb ge st ct el tr gn = 
    { chrom=chr; gbeg=gb; gend=ge; str=st; cat=ct;  exlist=el; trid=tr; gnid=gn}

  let chrom t = t.chrom
  let gbeg t = t.gbeg
  let gend t = t.gend
  let str t = t.str	
  let cat t = t.cat
  let exlist t = t.exlist
  let trid t = t.trid
  let gnid t = t.gnid
end



module Gene = 
struct
  type t = {
	chrom: string;
	gbeg: int;
	gend: int;
	str: strand;
	trarr: Transcript.t array;
	exarr: Exon.t array;
	gnid: string;
  }

  let create chr gb ge st ta ea gn = 
    { chrom=chr; gbeg=gb; gend=ge; str=st; trarr=ta; exarr= ea; gnid=gn}
	  
  let chrom g = g.chrom
  let gbeg g = g.gbeg
  let gend g = g.gend
  let str g = g.str	
  let trarr g = g.trarr
  let exarr g = g.exarr
  let gnid g = g.gnid
end



(* Exon projection. In a gene it is a maximal set of overlapping exons *)
module ExonProj =
struct
  type t = {
	chrom: string;
	gbeg: int;
	gend: int;
	str: strand;
	nbex: int;  (* number of exons if is made of *)
	exarr: Exon.t array;  (* The exons it is made of *)
	noingn: int;    (* number of the exon projection among all the exon projections of the gene (from 5') *)
	begingn: int;   (* begining of the exon projection in the virtual cDNA made by joining all the exon projections of the gene *)
	endingn: int;   (* end of the exon projection in the virtual cDNA made by joining all the exon projections of the gene *)
 }
	  
  let create chr gb ge st ne ea no bgn egn = 
    { chrom=chr; gbeg=gb; gend=ge; str=st; nbex=ne; exarr=ea; noingn=no; begingn=bgn; endingn=egn}
  
  let chrom e = e.chrom
  let gbeg e = e.gbeg
  let gend e = e.gend
  let str e = e.str
  let nbex e = e.nbex
  let exarr e = e.exarr
  let noingn e = e.noingn
  let begingn e = e.begingn 
  let endingn e = e.endingn
end 




(* Segmented (Exon) projection. The different exons of an exon projection define 
   elementary segments called Segmented (Exon) projections.
   The property of a segmented projection is that all its nt are included in the same
   number of transcripts.
*)
module SegProj = 
struct
  type t = {
	chrom: string;
	gbeg: int;
	gend: int;
	str: strand;
	score : float;
	exproj: ExonProj.t;
	coverage: int;   (* number of transcripts it belongs to *)
	noinexproj: int; (* number of the segmented projection in the exon projection from the 5' *)
	begingn: int; (* begining of the segmented projection in the virtual cdna of the gene *)
	endingn: int; (* end of the segmented projection in the virtual cdna of the gene *)
	
 }

  let create chr gb ge st sc ep cov no bgn egn = 
    { chrom=chr; gbeg=gb; gend=ge; str=st; score=sc; exproj=ep; coverage=cov; noinexproj=no; begingn=bgn; endingn=egn}
  
  let setcoverage cov s = {s with coverage=cov}
  let setscore sc s = {s with score=sc}

  let chrom s = s.chrom
  let gbeg s = s.gbeg
  let gend s = s.gend
  let str s = s.str
  let score s = s.score
  let exproj s = s.exproj
  let coverage s = s.coverage
  let noinexproj s = s.noinexproj
  let begingn s = s.begingn
  let endingn s = s.endingn
end



