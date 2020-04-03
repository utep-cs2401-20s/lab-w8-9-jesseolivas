import static org.junit.Assert.assertEquals;
import static org.junit.jupiter.api.Assertions.*;
import org.junit.Assert;
import org.junit.jupiter.api.Test;

public class AminoAcidLLTester {

    @Test
    // Tests isSorted method with an unsorted list
    public void isSortedT1(){
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("GUGGUAAUGAUAAUACUGUUUUUCUCGUAA");
        assertEquals(false, test.isSorted());
    }

    @Test
    // Tests isSorted method with a sorted list
    public void isSortedT2(){
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("GCGGCAUGUGACGAUGAAUUUGGGGGC");
        assertEquals(true, test.isSorted());
    }

    @Test
    // Tests sort method by passing through sort and then checks if it worked by passing through
    // isSorted method for easier testing
    public void sortT1(){
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("UUUUUUUUGUUACUUAUGAUACAGCAAUAA");
        test = AminoAcidLL.sort(test);
        assertEquals(true, test.isSorted());

    }

    @Test
    // Tests sort method by passing through sort and then checks if it worked by passing through
    // isSorted method for easier testing but with a one codon sequence
    public void sortT2(){
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("UUU");
        test = AminoAcidLL.sort(test);
        assertEquals(true, test.isSorted());

    }

    @Test
    // Tests createFromRNASequence to see if it works
    public void createFromRNASequenceT1(){
        String expected = "PQRVWST";
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("CCCCAGCGUGUGUGGAGCACGUAG");
        for (int i = 0; i < expected.length(); i++) {
            assertEquals(expected.charAt(i), test.aminoAcid);
            test = test.next;
        }
    }

    @Test
    // Tests createFromRNASequence with stop codon in the middle
    public void createFromRNASequenceT2(){
        String expected = "PQR";
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("CCCCAGCGUUAGGUGUGGAGCACGUAG");
        for (int i = 0; i < expected.length(); i++) {
            assertEquals(expected.charAt(i), test.aminoAcid);
            test = test.next;
        }
    }


}
