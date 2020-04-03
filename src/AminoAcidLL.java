class AminoAcidLL{
  char aminoAcid;
  String[] codons;
  int[] counts;
  AminoAcidLL next;

  AminoAcidLL(){

  }


  /********************************************************************************************/
  /* Creates a new node, with a given amino acid/codon 
   * pair and increments the codon counter for that codon.
   * NOTE: Does not check for repeats!! */
  AminoAcidLL(String inCodon){
    aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
    codons = AminoAcidResources.getCodonListForAminoAcid(aminoAcid);
    counts = new int[codons.length];

    // Counts how many codons made this specific amino acid
    incrementCodons(inCodon);
    next = null;
  }

  // Traverse through the codons and if they're the same,
  // increment at that index
  private void incrementCodons(String c){
    for (int i = 0; i < codons.length; i++){
      if(codons[i].equals(c)){
        counts[i]++;
      }
    }

  }
  /********************************************************************************************/
  /* Recursive method that increments the count for a specific codon:
   * If it should be at this node, increments it and stops, 
   * if not passes the task to the next node. 
   * If there is no next node, add a new node to the list that would contain the codon. 
   */
  private void addCodon(String inCodon){
    // if this node had this codon, increment count at that codon
    if(aminoAcid == AminoAcidResources.getAminoAcidFromCodon(inCodon)){
      incrementCodons(inCodon);
    }

    else {
      if (next != null) {
        // add a new node to the list that would contain the codon
        next.addCodon(inCodon);
      }
      else{
        next = new AminoAcidLL(inCodon);

      }
    }
  }


  /********************************************************************************************/
  /* Shortcut to find the total number of instances of this amino acid */
  private int totalCount(){
    int sum = 0;
    for(int i = 0; i < counts.length; i++){
      sum+= counts[i];
    }
    return sum;
  }

  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int totalDiff(AminoAcidLL inList){
    return Math.abs(totalCount() - inList.totalCount());
  }


  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int codonDiff(AminoAcidLL inList){
    int diff = 0;
    for(int i=0; i<codons.length; i++){
      diff += Math.abs(counts[i] - inList.counts[i]);
    }
    return diff;
  }

  /********************************************************************************************/
  /* Recursive method that finds the differences in **Amino Acid** counts. 
   * the list *must* be sorted to use this method */
  public int aminoAcidCompare(AminoAcidLL inList){
    // First sort the list
    inList = sort(inList);

    // Once both are at their end, you're done and return
    if(inList.next == null && next == null){
      return codonDiff(inList);
    }
    else if(next == null){
      return inList.totalCount() + aminoAcidCompare(inList);
    }
    else if(inList.next == null){
      return totalCount() + next.aminoAcidCompare(inList);
    }
    // Go to the next comparison
    return totalDiff(inList) + next.aminoAcidCompare(inList.next);
  }

  /********************************************************************************************/
  /* Same as above, but counts the codon usage differences
   * Must be sorted. */
  public int codonCompare(AminoAcidLL inList){
    // First sort the list
    inList = sort(inList);

    // Once both are at their end, you're done and return
    if(inList.next == null && next == null){
      return codonDiff(inList);
    }
    else if(next == null){
      return inList.totalCount() + codonCompare(inList);
    }
    else if(inList.next == null){
      return totalCount() + next.codonCompare(inList);
    }
    // Go to the next comparison
    return totalDiff(inList) + next.codonCompare(inList.next);
  }


  /********************************************************************************************/
  /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
  public char[] aminoAcidList(){
    // base, if next is null
    if(next == null){
      return new char[aminoAcid];
    }

    // recursion, when next is not null
    char[] a = next.aminoAcidList();
    char[] ret = new char[a.length + 1];
    ret[0] = aminoAcid;
    for (int i = 1; i < a.length; i++) {
      ret[i] = a[i];
    }

    // Return the array of chars holding the amino acid list
    return ret;
  }

  /********************************************************************************************/
  /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
  public int[] aminoAcidCounts(){
    // base, if next is null
    if(next == null){
      return new int[totalCount()];
    }

    // recursion, when next is not null
    int[] a = next.aminoAcidCounts();
    int[] ret = new int[a.length+1];
    ret[0] = totalCount();
    for(int i = 1; i < a.length; i++){
      ret[i] = a[i];
    }

    // Return the array of int holding the total counts
    return ret;
  }


  /********************************************************************************************/
  /* recursively determines if a linked list is sorted or not */
  public boolean isSorted(){
    // Ff reached the end of the list, sorted
    if(next == null){
      return true;
    }

    // If the next node in the list is greater by ascii value,
    // it's not sorted so return false
    else if(aminoAcid > next.aminoAcid){
      return false;
    }

    // If everything is sorted so far, iterate to the next node
    return next.isSorted();
  }


  /********************************************************************************************/
  /* Static method for generating a linked list from an RNA sequence */
  public static AminoAcidLL createFromRNASequence(String inSequence){
    AminoAcidLL list = new AminoAcidLL(inSequence.substring(0,3));
    while(inSequence.length() > 3 && AminoAcidResources.getAminoAcidFromCodon(inSequence.substring(0,3)) != '*'){
      inSequence = inSequence.substring(3);
      list.addCodon(inSequence.substring(0,3));
    }
    return list;
  }


  /********************************************************************************************/
  /* sorts a list by amino acid character*/
  public static AminoAcidLL sort(AminoAcidLL inList){
    // Returns automatically if list is empty or only has one
    if(inList == null || inList.next == null){
      return inList;
    }

    // Sort list by selection sort
    char temp;
    for( AminoAcidLL i = inList; i.next != null; i = i.next){
      for (AminoAcidLL j = i.next; j != null; j = j.next){
        if(i.aminoAcid > j.aminoAcid){
          temp = i.aminoAcid;
          i.aminoAcid = j.aminoAcid;
          j.aminoAcid = temp;

        }
      }
    }
    return inList;
  }
}