
(* main.ml: Given a bedpe file of generalized exon to exon junctions with two 24 mers around the "splice sites",
   searches for duplicated sequences in these two 24 mers and outputs the longest one(s).
*)

open Common
open Config
open Feature
open Input
open Output
open Printf



(* duplseq is a duplicated sequence found in the 24mer of the don sequence donseq
   "close to" don site, and also found anywhere in the 24mer of the acc sequence.
   duplseq is a triplet (lg,donst,accst) where:
   - lg is the length of the duplicated sequence, 
   - donst is the start of the occurence in donseq,
   - accst is the staty of the occurence in accseq.
   This function returns a duplicated sequence extended from duplseq on the left 
   and on the right, ie which includes duplseq but that could be longer.
   - : string -> string -> 'a * int * int -> int * int * int = <fun>
   to test it:
   ----------
   let donseq = "ATCTGAGCCTGGGTACGTGCGCTC" and accseq = "TTTTGTTTTTAGCCTGGTATAAGC" and duplseq = (5,7,12);; 
   must give (5,10) on the left and (7,5,10) on both sides
   let donseq = "ATCTGAGCCTGGGTACGTGCGCTC" and accseq = "TTTTGATCTGAGCCTGGTATAAGC" and duplseq = (5,7,12);;
   must give (0,5) on the left
   let donseq = "ATCTGAGCCTGGGTACGTGCGCTC" and accseq = "TTTTGATCTGAGCCTGGGTACGAC" and duplseq = (5,7,12);;
   must give (16,21) on the right and (17,0,5) on both sides
*)
let extend_max_both_sides donseq accseq duplseq =
  let (lginit,donst,accst) = duplseq in
    (* first try to extend on the left *)
  let idonl = ref (donst-1) and iaccl = ref (accst-1) in
    while ((!idonl>=0) && (!iaccl>=0) && (((String.get donseq (!idonl))=(String.get accseq (!iaccl))))) do
      decr idonl;
      decr iaccl;
    done;  (* !idonl+1 and !iaccl+1 contain the values that we want *)
    (* then try to extend on the right *)
    let idonr = ref (donst+5) and iaccr = ref (accst+5) in
    while ((!idonr<=23) && (!iaccr<=23) && (((String.get donseq (!idonr))=(String.get accseq (!iaccr))))) do
      incr idonr;
      incr iaccr;
    done;  (* !idonr-1 and !iaccr-1 contain the values that we want *)
      ((!idonr-1-(!idonl+1)+1),!idonl+1,!iaccl+1);;



(* 
   function to compare two extended duplicated sequences present in donor and
   acceptor 24mer sequences.
   We want the longest one to be first, and then if the length is the same we
   want the one that starts first in the donor sequence.
   In case of equality at this level we want to put first the ones that starts 
   first in the acceptor sequence.
   val compare_duplseq : 'a * 'b * 'c -> 'a * 'b * 'c -> int = <fun>
*)
let compare_duplseq (lg1,std1,sta1) (lg2,std2,sta2) =
  if (lg1>lg2) then
    -1
  else
    begin
      if (lg1<lg2) then
	1
      else
	begin
	  if (std1<std2) then
	    -1
	  else
	    begin
	      if (std1>std2) then
		1
	      else
		begin
		  if (sta1<sta2) then
		    -1
		  else
		    begin
		      if (sta1>sta2) then
			1
		      else
			0
		    end
		end
	    end
	end
    end
    

(* This is the most important function used by the main function find_dupl_seq_in_introns().
   ***************************************************************************************
   donorseq is a 24mer around the donor site represented as a string (here no check for
   correctness of the alphabet) = 12bp in the 5' exon then 2 bp for donor (usually GT) 
   and finally 10bp outside the 5' exon.
   accseq is a 24mer around the acc site represented as a string (here no check for
   correctness of the alphabet) = 10 bp outside the 3' exon then 2bp for acceptor (usually AG)
   and finally 12bp in the 3' exon.
   Reports for each 5-mer "close to" the donor site (2bp) ("close to" = either overlapping 
   or adjacent to the 2bp donor site(2bp)), each occurence of this 5mer in the acc sequence 
   extended as much as possible on both sides.
   - : string -> string -> (int * int * int) list = <fun>

   Note: 5 can be taken out as a parameter afterwards
   Note2: be careful: we need to remove the redundancy here since we can have duplicates
   in case of extension from two overlapping donor anchors
*)
let find_don_anch_in_acc donseq accseq =
  (* 
     Define the donor duplicated sequence filter =
     the donor duplicated sequence must either be adjacent or overlapping the 2bp junction/region.
     Means that the donor duplicated sequence can start at one of the following positions: 6, 7, 8, 9, 10, 11, 12, 13 (0-based).
     Note:
     ****
     - 13 and 14 added after Jon's email of 08/06/2008:
     [CCTCAAGAAGG/gctacat...  ...gctcaagaaag/ACACTGATT], which has CTCAAGAA duplicated, must be detected
     - 14 removed and 6 added after Jon's email of 08/07/2008 
     HE defines the 2bp region as the 2bp of the junction, not the splice site.
     Used for the final output of the present function.
  *)
  let anchors_in_donor = List.map (fun pos -> (5,pos,-1)) [6;7;8;9;10;11;12;13] in 
  
  (* 
     Define the acceptor duplicated sequence filter = 
     the acceptor duplicated sequence must either be adjacent or overlapping the 2bp junction/region.
     Means that the interval defined by the duplicated sequence in the acceptor 24mer
     must include one of the following positions: 10, 11, 12 or 13 (0-based).
     Note:
     ****
     - 9 and 10 added after Jon's email of 08/06/08:	
     [CCTCAAGAAGG/gctacat...  ...gctcaagaaag/ACACTGATT], which has CTCAAGAA duplicated, must be detected 
     - 9 removed and 13 added after Jon's email of 08/07/2008:
     he defines the 2bp region as the 2bp of the junction, not the splice site.
  *)
  let accfilter (l,sd,sa) = ((Common.in_interval 10 (sa,sa+l-1))||(Common.in_interval 11 (sa,sa+l-1))||(Common.in_interval 12 (sa,sa+l-1))||(Common.in_interval 13 (sa,sa+l-1))) in
    

  (*
    Define the acceptor-donor filter = we cannot have any bp of the duplicated sequence that are not exonic on both sides.
    Being outside the exon for a donor position is being >= 12;
    Being outside the exon for an acceptor position is being <= 11. 
    Used for the final output of the present function.
  *)
  let donaccfilter (l,sd,sa) =
    let id = ref sd and ia = ref sa and ok= ref true in
    while (!id<=sd+l-1) do
      if((!id>=12)&&(!ia<=11)) then
	ok:= false;
      incr id;
      incr ia;
    done;
      !ok
  in
    
  (* Since we do not know how it will extend, we look at an acc duplicated sequence starting
     anywhere in the acc 24mer 
  *)
  let all5mersinacc = List.map (fun pos -> (5,pos)) [0;1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19] in
  let anchors_in_donor_with_occinacc = ref [] and lextended_duplseq = ref [] in
  
  (* For each anchor in donor site, report if there is a substring of the acc that matches *)
  let u = List.iter
    (
      fun (lgd,std,_) -> 
	List.iter 
	  (
	    fun (lga,sta) -> 
	      if (((String.compare (String.sub donseq std lgd) (String.sub accseq sta lga))) = 0) then 
		anchors_in_donor_with_occinacc:=(Common.insert_right_place_without_redond (lgd,std,sta) (!anchors_in_donor_with_occinacc) compare_duplseq)
	  )
	  all5mersinacc
    )
    anchors_in_donor in
    
  
    (* then for each exact match try to extend on both sides, first on the left, then on the right 
       and then apply the acc filter regarding position 
       and also the donacc filter regarding the rule of the intron = one bp from the DS cannot be in two introns
    *)
    List.filter donaccfilter (List.filter accfilter (remove_redund compare_duplseq (List.sort compare_duplseq (List.map (extend_max_both_sides donseq accseq) (!anchors_in_donor_with_occinacc)))));;


(*    
      not done yet 
      let select_longest_dupl_seq donseq accseq =
*) 



(* 
   Given a bedpe file of generalized exon to exon junctions with two 24 mers around the splice sites,
   outputs all the duplicated sequences included in these 24 mers which
   don and acc parts are "close to" the 2bp splice site on each side.
   "close to" = either overlapping or adjacent to the 2bp junction
*)
let find_dupl_seq_in_introns () =
  (* read the command line *)
  read_commandline ();
  
  (* read the bedpe file and store the info in a list of features *)
  let flist = make_feat_list context.file in
 
  (* apply the algorithm *)
  let flist2 = List.map 
    (fun feat -> 
      let donseq = List.nth (Feature.attlist feat) 0 and accseq = List.nth (Feature.attlist feat) 1 in
	    (* be careful: here we need to filter out the dupl sequences which acc part is not "close" to 2bp acc 
	       "close to" = either overlapping or adjacent to the 2bp junction.
	       Also we want to report the length and the start positions in each of the 24 mers, of the duplicated sequence
	    *)
	  let lduplseq = ref [] and llength = ref [] and lstartdon = ref [] and lstartacc = ref [] in
	  let u = List.iter 
	    (fun (lg,std,sta) -> 
		begin
		  lduplseq:=List.append (!lduplseq) [String.sub donseq std lg];
		  llength:=List.append (!llength) [string_of_int lg];
		  lstartdon:=List.append (!lstartdon) [string_of_int (std+1)];
		  lstartacc:=List.append (!lstartacc) [string_of_int (sta+1)];
		end
	    )
	    (find_don_anch_in_acc donseq accseq) 
	  in
	  let string_of_duplseq = Common.list_to_string (!lduplseq) and string_of_llength = Common.list_to_string (!llength) in
	  let string_of_lstartdon = Common.list_to_string (!lstartdon) and string_of_lstartacc = Common.list_to_string (!lstartacc) in
	    Feature.setattlist (List.append (Feature.attlist feat) [string_of_duplseq;string_of_llength;string_of_lstartdon;string_of_lstartacc]) feat
    )
    flist
  in

  let u = Common.print_log ((Filename.basename (Sys.argv.(0)))^(" did its work ! Now writing output file\n")) in

  (* write the output in context.outfile *)
  let o = try
      (open_out context.outfile) 
    with
      | Sys_error s -> open_out ((Sys.argv.(1))^("_withduplseq.bedpe"))
  in
  let u = print o flist2 in
  let u = flush o in
    ();;


find_dupl_seq_in_introns ();;
