package org.opencb.cellbase.core.variant.annotation;

import org.opencb.biodata.models.core.Exon;
import org.opencb.biodata.models.core.Gene;
import org.opencb.biodata.models.core.RegulatoryFeature;
import org.opencb.biodata.models.core.Transcript;
import org.opencb.biodata.models.variant.Variant;
import org.opencb.biodata.models.variant.avro.ConsequenceType;
import org.opencb.biodata.models.variant.avro.ProteinVariantAnnotation;
import org.opencb.cellbase.core.api.GenomeDBAdaptor;
import org.opencb.commons.datastore.core.Query;
import org.opencb.commons.datastore.core.QueryOptions;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

/**
 * Created by fjlopez on 15/03/17.
 */
public class ConsequenceTypeMNVCalculator extends ConsequenceTypeCalculator {
    private int variantStart;
    private int variantEnd;
    private GenomeDBAdaptor genomeDBAdaptor;

    private static final int BIG_VARIANT_SIZE_THRESHOLD = 50;


    public ConsequenceTypeMNVCalculator(GenomeDBAdaptor genomeDBAdaptor) {
        this.genomeDBAdaptor = genomeDBAdaptor;
    }

    public List<ConsequenceType> run(Variant inputVariant, List<Gene> geneList, List<RegulatoryFeature> regulatoryRegionList) {
        List<ConsequenceType> consequenceTypeList = new ArrayList<>();
        variant = inputVariant;
        variantEnd = variant.getStart() + variant.getReference().length() - 1;
        variantStart = variant.getStart();
        boolean isIntergenic = true;
        for (Gene currentGene : geneList) {
            gene = currentGene;
            for (Transcript currentTranscript : gene.getTranscripts()) {
                isIntergenic = isIntergenic && (variantEnd < currentTranscript.getStart() || variantStart > currentTranscript.getEnd());
                transcript = currentTranscript;
                consequenceType = new ConsequenceType();
                consequenceType.setGeneName(gene.getName());
                consequenceType.setEnsemblGeneId(gene.getId());
                consequenceType.setEnsemblTranscriptId(transcript.getId());
                consequenceType.setStrand(transcript.getStrand());
                consequenceType.setBiotype(transcript.getBiotype());
                consequenceType.setTranscriptAnnotationFlags(transcript.getAnnotationFlags() != null
                        ? new ArrayList<>(transcript.getAnnotationFlags()) : null);
                SoNames.clear();

                if (transcript.getStrand().equals("+")) {
                    // Deletion - whole transcript removed
                    if (variantStart <= transcript.getStart() && variantEnd >= transcript.getEnd()) {
                        SoNames.add(VariantAnnotationUtils.TRANSCRIPT_ABLATION);
                        consequenceType.setSequenceOntologyTerms(getSequenceOntologyTerms(SoNames));
                        consequenceTypeList.add(consequenceType);
                    } else if (regionsOverlap(transcript.getStart(), transcript.getEnd(), variantStart, variantEnd)) {
                        switch (transcript.getBiotype()) {
                            /**
                             * Coding biotypes
                             */
                            case VariantAnnotationUtils.NONSENSE_MEDIATED_DECAY:
                                SoNames.add(VariantAnnotationUtils.NMD_TRANSCRIPT_VARIANT);
                            case VariantAnnotationUtils.IG_C_GENE:
                            case VariantAnnotationUtils.IG_D_GENE:
                            case VariantAnnotationUtils.IG_J_GENE:
                            case VariantAnnotationUtils.IG_V_GENE:
                            case VariantAnnotationUtils.TR_C_GENE:  // TR_C_gene
                            case VariantAnnotationUtils.TR_D_GENE:  // TR_D_gene
                            case VariantAnnotationUtils.TR_J_GENE:  // TR_J_gene
                            case VariantAnnotationUtils.TR_V_GENE:  // TR_V_gene
                            case VariantAnnotationUtils.POLYMORPHIC_PSEUDOGENE:
                            case VariantAnnotationUtils.PROTEIN_CODING:    // protein_coding
                            case VariantAnnotationUtils.NON_STOP_DECAY:    // non_stop_decay
                            case VariantAnnotationUtils.TRANSLATED_PROCESSED_PSEUDOGENE:
                            case VariantAnnotationUtils.TRANSLATED_UNPROCESSED_PSEUDOGENE:    // translated_unprocessed_pseudogene
                            case VariantAnnotationUtils.LRG_GENE:    // LRG_gene
                                solveCodingPositiveTranscript();
                                consequenceType.setSequenceOntologyTerms(getSequenceOntologyTerms(SoNames));
                                consequenceTypeList.add(consequenceType);
                                break;
                            /**
                             * Non-coding biotypes
                             */
                            default:
                                solveNonCodingPositiveTranscript();
                                consequenceType.setSequenceOntologyTerms(getSequenceOntologyTerms(SoNames));
                                consequenceTypeList.add(consequenceType);
                                break;
                        }
                    } else {
                        solveTranscriptFlankingRegions(VariantAnnotationUtils.UPSTREAM_GENE_VARIANT,
                                VariantAnnotationUtils.DOWNSTREAM_GENE_VARIANT);
                        if (SoNames.size() > 0) { // Variant does not overlap gene region, just may have upstream/downstream annotations
                            consequenceType.setSequenceOntologyTerms(getSequenceOntologyTerms(SoNames));
                            consequenceTypeList.add(consequenceType);
                        }
                    }
                } else {
                    if (variantStart <= transcript.getStart() && variantEnd >= transcript.getEnd()) { // Deletion - whole trans. removed
                        SoNames.add(VariantAnnotationUtils.TRANSCRIPT_ABLATION);
                        consequenceType.setSequenceOntologyTerms(getSequenceOntologyTerms(SoNames));
                        consequenceTypeList.add(consequenceType);
                    } else if (regionsOverlap(transcript.getStart(), transcript.getEnd(), variantStart, variantEnd)) {
                        switch (transcript.getBiotype()) {
                            /**
                             * Coding biotypes
                             */
                            case VariantAnnotationUtils.NONSENSE_MEDIATED_DECAY:
                                SoNames.add(VariantAnnotationUtils.NMD_TRANSCRIPT_VARIANT);
                            case VariantAnnotationUtils.IG_C_GENE:
                            case VariantAnnotationUtils.IG_D_GENE:
                            case VariantAnnotationUtils.IG_J_GENE:
                            case VariantAnnotationUtils.IG_V_GENE:
                            case VariantAnnotationUtils.TR_C_GENE:  // TR_C_gene
                            case VariantAnnotationUtils.TR_D_GENE:  // TR_D_gene
                            case VariantAnnotationUtils.TR_J_GENE:  // TR_J_gene
                            case VariantAnnotationUtils.TR_V_GENE:  // TR_V_gene
                            case VariantAnnotationUtils.POLYMORPHIC_PSEUDOGENE:
                            case VariantAnnotationUtils.PROTEIN_CODING:    // protein_coding
                            case VariantAnnotationUtils.NON_STOP_DECAY:    // non_stop_decay
                            case VariantAnnotationUtils.TRANSLATED_PROCESSED_PSEUDOGENE:
                            case VariantAnnotationUtils.TRANSLATED_UNPROCESSED_PSEUDOGENE:  // translated_unprocessed_pseudogene
                            case VariantAnnotationUtils.LRG_GENE:    // LRG_gene
                                solveCodingNegativeTranscript();
                                consequenceType.setSequenceOntologyTerms(getSequenceOntologyTerms(SoNames));
                                consequenceTypeList.add(consequenceType);
                                break;
                            /**
                             * Non-coding biotypes
                             */
                            default:
                                solveNonCodingNegativeTranscript();
                                consequenceType.setSequenceOntologyTerms(getSequenceOntologyTerms(SoNames));
                                consequenceTypeList.add(consequenceType);
                                break;
                        }
                    } else {
                        solveTranscriptFlankingRegions(VariantAnnotationUtils.DOWNSTREAM_GENE_VARIANT,
                                VariantAnnotationUtils.UPSTREAM_GENE_VARIANT);
                        if (SoNames.size() > 0) { // Variant does not overlap gene region, just has upstream/downstream annotations
                            consequenceType.setSequenceOntologyTerms(getSequenceOntologyTerms(SoNames));
                            consequenceTypeList.add(consequenceType);
                        }
                    }
                }
            }
        }

        if (consequenceTypeList.size() == 0 && isIntergenic) {
            HashSet<String> intergenicName = new HashSet<>();
            intergenicName.add(VariantAnnotationUtils.INTERGENIC_VARIANT);
            ConsequenceType consequenceType = new ConsequenceType();
            consequenceType.setSequenceOntologyTerms(getSequenceOntologyTerms(intergenicName));
            consequenceTypeList.add(consequenceType);
        }

        solveRegulatoryRegions(regulatoryRegionList, consequenceTypeList);
        return consequenceTypeList;
    }

    private void solveNonCodingNegativeTranscript() {
        Exon exon = transcript.getExons().get(0);
        boolean variantAhead = true; // we need a first iteration within the while to ensure junction is solved in case needed
        int cdnaExonEnd = (exon.getEnd() - exon.getStart() + 1);
        int cdnaVariantStart = -1;
        int cdnaVariantEnd = -1;
        int firstCdsPhase = -1;
        boolean[] junctionSolution = {false, false};
        boolean splicing = false;
        Integer variantStartExonNumber = null;
        Integer variantEndExonNumber = null;

        if (firstCdsPhase == -1 && transcript.getGenomicCodingEnd() >= exon.getStart()) {
            firstCdsPhase = exon.getPhase();
        }
        if (variantEnd <= exon.getEnd()) {
            if (variantEnd >= exon.getStart()) {  // Variant end within the exon
                cdnaVariantStart = cdnaExonEnd - (variantEnd - exon.getStart());
                consequenceType.setCdnaPosition(cdnaVariantStart);
                variantStartExonNumber = exon.getExonNumber();
                if (variantStart >= exon.getStart()) {  // Both variant start and variant end within the exon  ----||||S|||||E||||----
                    cdnaVariantEnd = cdnaExonEnd - (variantStart - exon.getStart());
                }
            }
        } else if (variantStart >= exon.getStart()) {
            // We do not contemplate that variant end can be located before this exon since this is the first exon
            cdnaVariantEnd = cdnaExonEnd - (variantEnd - exon.getStart());
            variantEndExonNumber = exon.getExonNumber();
        } // Variant includes the whole exon. Variant end is located before the exon, variant start is located after the exon

        int exonCounter = 1;
        // This is not a do-while since we cannot call solveJunction until
        while (exonCounter < transcript.getExons().size() && variantAhead) {
            int prevSpliceSite = exon.getStart() - 1;
            exon = transcript.getExons().get(exonCounter);          // next exon has been loaded
            // Set firsCdsPhase only when the first coding exon is reached
            if (firstCdsPhase == -1 && transcript.getGenomicCodingEnd() >= exon.getStart()) {
                firstCdsPhase = exon.getPhase();
            }
            solveJunction(exon.getEnd() + 1, prevSpliceSite, VariantAnnotationUtils.SPLICE_ACCEPTOR_VARIANT,
                    VariantAnnotationUtils.SPLICE_DONOR_VARIANT, junctionSolution);
            splicing = (splicing || junctionSolution[0]);

            if (variantEnd <= exon.getEnd()) {
                cdnaExonEnd += (exon.getEnd() - exon.getStart() + 1);
                if (variantEnd >= exon.getStart()) {  // Variant end within the exon
                    cdnaVariantStart = cdnaExonEnd - (variantEnd - exon.getStart());
                    consequenceType.setCdnaPosition(cdnaVariantStart);
                    variantStartExonNumber = exon.getExonNumber();
                    if (variantStart >= exon.getStart()) {  // Both variant start and variant end within the exon  ----||||S|||||E||||----
                        cdnaVariantEnd = cdnaExonEnd - (variantStart - exon.getStart());
                    }
                }
            } else if (variantStart >= exon.getStart()) {
                if (variantStart <= exon.getEnd()) {  // Only variant start within the exon  ----||||||||||E||||----
                    cdnaExonEnd += (exon.getEnd() - exon.getStart() + 1);
                    cdnaVariantEnd = cdnaExonEnd - (variantStart - exon.getStart());
                    variantEndExonNumber = exon.getExonNumber();
                } else {  // Variant does not include this exon, variant is located before this exon
                    variantAhead = false;
                }
            } else {  // Variant includes the whole exon. Variant start is located before the exon, variant end is located after the exon
                cdnaExonEnd += (exon.getEnd() - exon.getStart() + 1);
            }
            exonCounter++;
        }
        consequenceType.setExonNumber(variantStartExonNumber != null ? variantStartExonNumber : variantEndExonNumber);
        solveMiRNA(cdnaVariantStart, cdnaVariantEnd, junctionSolution[1]);
    }

    private void solveCodingNegativeTranscript() {
        Exon exon = transcript.getExons().get(0);
        String transcriptSequence = exon.getSequence();
        boolean variantAhead = true; // we need a first iteration within the while to ensure junction is solved in case needed
        int cdnaExonEnd = (exon.getEnd() - exon.getStart() + 1);
        int cdnaVariantStart = -1;
        int cdnaVariantEnd = -1;
        int firstCdsPhase = -1;
        boolean[] junctionSolution = {false, false};
        boolean splicing = false;
        Integer variantStartExonNumber = null;
        Integer variantEndExonNumber = null;

        if (firstCdsPhase == -1 && transcript.getGenomicCodingEnd() >= exon.getStart()) {
            firstCdsPhase = exon.getPhase();
        }
        if (variantEnd <= exon.getEnd()) {
            if (variantEnd >= exon.getStart()) {  // Variant end within the exon
                cdnaVariantStart = cdnaExonEnd - (variantEnd - exon.getStart());
                consequenceType.setCdnaPosition(cdnaVariantStart);
                variantStartExonNumber = exon.getExonNumber();
                if (variantStart >= exon.getStart()) {  // Both variant start and variant end within the exon  ----||||S|||||E||||----
                    cdnaVariantEnd = cdnaExonEnd - (variantStart - exon.getStart());
                }
            }
        } else if (variantStart >= exon.getStart()) {
            // We do not contemplate that variant end can be located before this exon since this is the first exon
            cdnaVariantEnd = cdnaExonEnd - (variantEnd - exon.getStart());
            variantEndExonNumber = exon.getExonNumber();
        } // Variant includes the whole exon. Variant end is located before the exon, variant start is located after the exon

        int exonCounter = 1;
        // This is not a do-while since we cannot call solveJunction  until
        while (exonCounter < transcript.getExons().size() && variantAhead) {
            int prevSpliceSite = exon.getStart() - 1;
            exon = transcript.getExons().get(exonCounter);          // next exon has been loaded
            transcriptSequence = exon.getSequence() + transcriptSequence;
            // Set firsCdsPhase only when the first coding exon is reached
            if (firstCdsPhase == -1 && transcript.getGenomicCodingEnd() >= exon.getStart()) {
                firstCdsPhase = exon.getPhase();
            }
            solveJunction(exon.getEnd() + 1, prevSpliceSite, VariantAnnotationUtils.SPLICE_ACCEPTOR_VARIANT,
                    VariantAnnotationUtils.SPLICE_DONOR_VARIANT, junctionSolution);
            splicing = (splicing || junctionSolution[0]);

            if (variantEnd <= exon.getEnd()) {
                cdnaExonEnd += (exon.getEnd() - exon.getStart() + 1);
                if (variantEnd >= exon.getStart()) {  // Variant end within the exon
                    cdnaVariantStart = cdnaExonEnd - (variantEnd - exon.getStart());
                    consequenceType.setCdnaPosition(cdnaVariantStart);
                    variantStartExonNumber = exon.getExonNumber();
                    if (variantStart >= exon.getStart()) {  // Both variant start and variant end within the exon  ----||||S|||||E||||----
                        cdnaVariantEnd = cdnaExonEnd - (variantStart - exon.getStart());
                    }
                }
            } else if (variantStart >= exon.getStart()) {
                if (variantStart <= exon.getEnd()) {  // Only variant start within the exon  ----||||||||||E||||----
                    cdnaExonEnd += (exon.getEnd() - exon.getStart() + 1);
                    cdnaVariantEnd = cdnaExonEnd - (variantStart - exon.getStart());
                    variantEndExonNumber = exon.getExonNumber();
                } else {  // Variant does not include this exon, variant is located before this exon
                    variantAhead = false;
                }
            } else {  // Variant includes the whole exon. Variant start is located before the exon, variant end is located after the exon
                cdnaExonEnd += (exon.getEnd() - exon.getStart() + 1);
            }
            exonCounter++;
        }

        consequenceType.setExonNumber(variantStartExonNumber != null ? variantStartExonNumber : variantEndExonNumber);
        // Is not intron variant (both ends fall within the same intron)
        if (!junctionSolution[1]) {
            solveExonVariantInNegativeTranscript(splicing, transcriptSequence, cdnaVariantStart, cdnaVariantEnd,
                    firstCdsPhase);
        }
    }

    private void solveExonVariantInNegativeTranscript(boolean splicing, String transcriptSequence,
                                                      int cdnaVariantStart, int cdnaVariantEnd, int firstCdsPhase) {
        if (variantEnd > transcript.getGenomicCodingEnd()) {
            if (transcript.getEnd() > transcript.getGenomicCodingEnd() || transcript.unconfirmedStart()) { // Check transcript has 3 UTR
                SoNames.add(VariantAnnotationUtils.FIVE_PRIME_UTR_VARIANT);
            }
            if (variantStart <= transcript.getGenomicCodingEnd()) {
                SoNames.add(VariantAnnotationUtils.CODING_SEQUENCE_VARIANT);
                // cdnaCodingStart < 1 if cds_start_NF and phase!=0
                if (transcript.getCdnaCodingStart() > 0 || !transcript.unconfirmedStart()) {
                    SoNames.add(VariantAnnotationUtils.START_LOST);
                }
                if (variantStart < (transcript.getGenomicCodingStart() + 3)) {
                    SoNames.add(VariantAnnotationUtils.STOP_LOST);
                    if (variantStart < transcript.getGenomicCodingStart()) {
                        if (transcript.getStart() < transcript.getGenomicCodingStart()
                                || transcript.unconfirmedEnd()) { // Check transcript has 3 UTR)
                            SoNames.add(VariantAnnotationUtils.THREE_PRIME_UTR_VARIANT);
                        }
                    }
                }
            }
        } else if (variantEnd >= transcript.getGenomicCodingStart()) {
            // Need to define a local cdnaCodingStart because may modified in two lines below
            int cdnaCodingStart = transcript.getCdnaCodingStart();
            if (cdnaVariantStart != -1) {  // cdnaVariantStart may be null if variantEnd falls in an intron
                if (transcript.unconfirmedStart()) {
                    cdnaCodingStart -= ((3 - firstCdsPhase) % 3);
                }
                int cdsVariantStart = cdnaVariantStart - cdnaCodingStart + 1;
                consequenceType.setCdsPosition(cdsVariantStart);
                // First place where protein variant annotation is added to the Consequence type,
                // must create the ProteinVariantAnnotation object
                ProteinVariantAnnotation proteinVariantAnnotation = new ProteinVariantAnnotation();
                proteinVariantAnnotation.setPosition(((cdsVariantStart - 1) / 3) + 1);
                consequenceType.setProteinVariantAnnotation(proteinVariantAnnotation);
            }
            if (variantStart >= transcript.getGenomicCodingStart()) {  // Variant start also within coding region
                solveCodingExonVariantInNegativeTranscript(splicing, transcriptSequence, cdnaCodingStart, cdnaVariantStart,
                        cdnaVariantEnd);
            } else {
                // Check transcript has 3 UTR)
                if (transcript.getStart() < transcript.getGenomicCodingStart() || transcript.unconfirmedEnd()) {
                    SoNames.add(VariantAnnotationUtils.THREE_PRIME_UTR_VARIANT);
                }
                SoNames.add(VariantAnnotationUtils.CODING_SEQUENCE_VARIANT);
                SoNames.add(VariantAnnotationUtils.STOP_LOST);
            }
            // Check transcript has 3 UTR)
        } else if (transcript.getStart() < transcript.getGenomicCodingStart() || transcript.unconfirmedEnd()) {
            SoNames.add(VariantAnnotationUtils.THREE_PRIME_UTR_VARIANT);
        }
    }

    private void solveCodingExonVariantInNegativeTranscript(boolean splicing, String transcriptSequence, int cdnaCodingStart,
                                                            int cdnaVariantStart, int cdnaVariantEnd) {
        Boolean codingAnnotationAdded = false;

        // cdnaVariantStart=null if variant is intronic. cdnaCodingStart<1 if cds_start_NF and phase!=0
        if (cdnaVariantStart != -1 && cdnaVariantStart < (cdnaCodingStart + 3) && (cdnaCodingStart > 0
                || !transcript.unconfirmedStart())) {
            SoNames.add(VariantAnnotationUtils.START_LOST);
            codingAnnotationAdded = true;
        }
        if (cdnaVariantEnd != -1) {
            int finalNtPhase = (transcript.getCdnaCodingEnd() - cdnaCodingStart) % 3;
            Boolean stopToSolve = true;
            // just checks cdnaVariantStart!=null because no splicing means cdnaVariantEnd is also != null
            if (!splicing && cdnaVariantStart != -1) {
                codingAnnotationAdded = true;
                if (Math.abs(variant.getReference().length() - variant.getAlternate().length()) % 3 == 0) {
                    SoNames.add(VariantAnnotationUtils.INFRAME_VARIANT);
                } else {
                    SoNames.add(VariantAnnotationUtils.FRAMESHIFT_VARIANT);
                }
                stopToSolve = false;  // Stop codon annotation will be solved in the line below.
                solveStopCodonNegativeMNV(transcriptSequence, cdnaCodingStart, cdnaVariantStart, cdnaVariantEnd);
            }
            if (cdnaVariantEnd >= (transcript.getCdnaCodingEnd() - finalNtPhase)) {
                if (finalNtPhase != 2) {
                    SoNames.add(VariantAnnotationUtils.INCOMPLETE_TERMINAL_CODON_VARIANT);
                } else if (stopToSolve) {  // Only if stop codon annotation was not already solved in the if block above
                    SoNames.add(VariantAnnotationUtils.STOP_LOST);
                }
            }
        }
        if (!codingAnnotationAdded) {
            SoNames.add(VariantAnnotationUtils.CODING_SEQUENCE_VARIANT);
        }
    }

    private void solveStopCodonNegativeMNV(String transcriptSequence, int cdnaCodingStart,
                                           int cdnaVariantStart, int cdnaVariantEnd) {
        Integer variantPhaseShift1 = (cdnaVariantStart - cdnaCodingStart) % 3;
        Integer variantPhaseShift2 = (cdnaVariantEnd - cdnaCodingStart) % 3;
        int modifiedCodon1Start = cdnaVariantStart - variantPhaseShift1;
        int modifiedCodon2Start = cdnaVariantEnd - variantPhaseShift2;
        if (modifiedCodon1Start > 0 && (modifiedCodon2Start + 2) <= transcriptSequence.length()) {
            /**
             * Calculate reference codon
             */
            String reverseCodon1 = new StringBuilder(transcriptSequence.substring(transcriptSequence.length() - modifiedCodon1Start - 2,
                    // Rigth limit of the substring sums +1 because substring does not include that position
                    transcriptSequence.length() - modifiedCodon1Start + 1)).reverse().toString();
            String reverseCodon2 = new StringBuilder(transcriptSequence.substring(transcriptSequence.length() - modifiedCodon2Start - 2,
                    // Rigth limit of the substring sums +1 because substring does not include that position
                    transcriptSequence.length() - modifiedCodon2Start + 1)).reverse().toString();
            String reverseTranscriptSequence = new StringBuilder(
                    transcriptSequence.substring(((transcriptSequence.length() - cdnaVariantEnd) > 2)
                                    ? (transcriptSequence.length() - cdnaVariantEnd - 3)
                                    : 0,  // Be careful reaching the end of the transcript sequence
                            // Rigth limit of the substring -2 because substring does not include that position
                            transcriptSequence.length() - cdnaVariantEnd)).reverse().toString();
            char[] referenceCodon1Array = reverseCodon1.toCharArray();
            referenceCodon1Array[0] = VariantAnnotationUtils.COMPLEMENTARY_NT.get(referenceCodon1Array[0]);
            referenceCodon1Array[1] = VariantAnnotationUtils.COMPLEMENTARY_NT.get(referenceCodon1Array[1]);
            referenceCodon1Array[2] = VariantAnnotationUtils.COMPLEMENTARY_NT.get(referenceCodon1Array[2]);
            String referenceCodon1 = String.valueOf(referenceCodon1Array);
            char[] referenceCodon2Array = reverseCodon2.toCharArray();
            referenceCodon2Array[0] = VariantAnnotationUtils.COMPLEMENTARY_NT.get(referenceCodon2Array[0]);
            referenceCodon2Array[1] = VariantAnnotationUtils.COMPLEMENTARY_NT.get(referenceCodon2Array[1]);
            referenceCodon2Array[2] = VariantAnnotationUtils.COMPLEMENTARY_NT.get(referenceCodon2Array[2]);
            String referenceCodon2 = String.valueOf(referenceCodon2Array);


            // Char array to contain the upper/lower-case formatted string for the codon change, e.g. aGT/ATG
            char[] formattedReferenceCodon1Array = String.valueOf(referenceCodon1Array).toLowerCase().toCharArray();


            /**
             * Calculate codons generated by the insertion
             */
            int i = 0;
            int reverseTranscriptSequencePosition = 0;
            char[] altArray = (new StringBuilder(variant.getAlternate()).reverse().toString()).toCharArray();
            int modifiedCodonPosition;
            int modifiedCodonPositionStart = variantPhaseShift1;
            char[] formattedModifiedCodonArray = referenceCodon1.toLowerCase().toCharArray();
            char[] modifiedCodonArray = referenceCodon1Array.clone();
            boolean useMitochondrialCode = variant.getChromosome().equals("MT");
            boolean firstCodon = true;

            do {
                for (modifiedCodonPosition = modifiedCodonPositionStart;
                    // Paste alternative nt in the corresponding codon position
                     (modifiedCodonPosition < 3 && i < variant.getAlternate().length()); modifiedCodonPosition++) {
                    modifiedCodonArray[modifiedCodonPosition] = VariantAnnotationUtils.COMPLEMENTARY_NT.get(altArray[i]);

                    // Edit modified nt to make it upper-case in the formatted strings
                    formattedReferenceCodon1Array[modifiedCodonPosition]
                            = Character.toUpperCase(formattedReferenceCodon1Array[modifiedCodonPosition]);
                    formattedModifiedCodonArray[modifiedCodonPosition]
                            = Character.toUpperCase(modifiedCodonArray[modifiedCodonPosition]);

                    i++;
                }
                for (; modifiedCodonPosition < 3; modifiedCodonPosition++) {  // Concatenate reference codon nts after alternative nts
                    if (reverseTranscriptSequencePosition >= reverseTranscriptSequence.length()) {
                        int genomicCoordinate = transcript.getStart()
                                - (reverseTranscriptSequencePosition - reverseTranscriptSequence.length() + 1);
                        Query query = new Query(GenomeDBAdaptor.QueryParams.REGION.key(), variant.getChromosome()
                                + ":" + genomicCoordinate
                                + "-" + (genomicCoordinate + 1));
                        modifiedCodonArray[modifiedCodonPosition] = VariantAnnotationUtils.COMPLEMENTARY_NT
                                .get(genomeDBAdaptor.getGenomicSequence(query, new QueryOptions())
                                        .getResult().get(0).getSequence().charAt(0));
                    } else {
                        modifiedCodonArray[modifiedCodonPosition] = VariantAnnotationUtils.COMPLEMENTARY_NT.get(
                                reverseTranscriptSequence.charAt(reverseTranscriptSequencePosition));
                    }
                    reverseTranscriptSequencePosition++;

                    // Edit modified nt to make it upper-case in the formatted strings
                    formattedReferenceCodon1Array[modifiedCodonPosition]
                            = Character.toUpperCase(formattedReferenceCodon1Array[modifiedCodonPosition]);
                    formattedModifiedCodonArray[modifiedCodonPosition]
                            = Character.toUpperCase(modifiedCodonArray[modifiedCodonPosition]);
                }

                // Set codon str, protein ref and protein alt ONLY for the first codon mofified by the insertion
                if (firstCodon) {
                    firstCodon = false;
                    // Only the exact codon where the deletion starts is set
                    consequenceType.setCodon(String.valueOf(formattedReferenceCodon1Array) + "/"
                            + String.valueOf(formattedModifiedCodonArray));
                    // Assumes proteinVariantAnnotation attribute is already initialized
                    consequenceType
                            .getProteinVariantAnnotation()
                            .setReference(VariantAnnotationUtils.getAminoacid(useMitochondrialCode, referenceCodon1));
                    consequenceType
                            .getProteinVariantAnnotation()
                            .setAlternate(VariantAnnotationUtils.getAminoacid(useMitochondrialCode,
                                    String.valueOf(modifiedCodonArray)));
                }

                decideStopCodonModificationAnnotation(SoNames,
                        VariantAnnotationUtils.isStopCodon(useMitochondrialCode, referenceCodon2)
                                ? referenceCodon2 : referenceCodon1, String.valueOf(modifiedCodonArray),
                        useMitochondrialCode);

                modifiedCodonPositionStart = 0;  // Reset the position where the next modified codon must be started to be filled
            } while (i < variant.getAlternate().length());  // All posible new codons generated by the inserted sequence must be checked

        }
    }

    private void solveTranscriptFlankingRegions(String leftRegionTag, String rightRegionTag) {
        // Variant overlaps with -5kb region
        if (regionsOverlap(transcript.getStart() - 5000, transcript.getStart() - 1, variantStart, variantEnd)) {
            // Variant overlaps with -2kb region
            if (regionsOverlap(transcript.getStart() - 2000, transcript.getStart() - 1, variantStart, variantEnd)) {
                SoNames.add("2KB_" + leftRegionTag);
            } else {
                SoNames.add(leftRegionTag);
            }
        }
        // Variant overlaps with +5kb region
        if (regionsOverlap(transcript.getEnd() + 1, transcript.getEnd() + 5000, variantStart, variantEnd)) {
            // Variant overlaps with +2kb region
            if (regionsOverlap(transcript.getEnd() + 1, transcript.getEnd() + 2000, variantStart, variantEnd)) {
                SoNames.add("2KB_" + rightRegionTag);
            } else {
                SoNames.add(rightRegionTag);
            }
        }
    }

    private void solveCodingPositiveTranscript() {
        Exon exon = transcript.getExons().get(0);
        String transcriptSequence = exon.getSequence();
        boolean variantAhead = true; // we need a first iteration within the while to ensure junction is solved in case needed
        int cdnaExonEnd = (exon.getEnd() - exon.getStart() + 1);
        int cdnaVariantStart = -1;
        int cdnaVariantEnd = -1;
        int firstCdsPhase = -1;
        boolean[] junctionSolution = {false, false};
        boolean splicing = false;
        Integer variantStartExonNumber = null;
        Integer variantEndExonNumber = null;

        if (transcript.getGenomicCodingStart() <= exon.getEnd()) {
            firstCdsPhase = exon.getPhase();
        }
        if (variantStart >= exon.getStart()) {
            if (variantStart <= exon.getEnd()) {  // Variant start within the exon
                cdnaVariantStart = cdnaExonEnd - (exon.getEnd() - variantStart);
                consequenceType.setCdnaPosition(cdnaVariantStart);
                variantStartExonNumber = exon.getExonNumber();
                if (variantEnd <= exon.getEnd()) {  // Both variant start and variant end within the exon  ----||||S|||||E||||----
                    cdnaVariantEnd = cdnaExonEnd - (exon.getEnd() - variantEnd);
                }
            }
        } else if (variantEnd <= exon.getEnd()) {
            // We do not contemplate that variant end can be located before this exon since this is the first exon
            cdnaVariantEnd = cdnaExonEnd - (exon.getEnd() - variantEnd);
            variantEndExonNumber = exon.getExonNumber();
        } // Variant includes the whole exon. Variant start is located before the exon, variant end is located after the exon


        int exonCounter = 1;
        // This is not a do-while since we cannot call solveJunction  until
        while (exonCounter < transcript.getExons().size() && variantAhead) {
            int prevSpliceSite = exon.getEnd() + 1;
            exon = transcript.getExons().get(exonCounter);          // next exon has been loaded
            transcriptSequence = transcriptSequence + exon.getSequence();
            // Set firsCdsPhase only when the first coding exon is reached
            if (firstCdsPhase == -1 && transcript.getGenomicCodingStart() <= exon.getEnd()) {
                firstCdsPhase = exon.getPhase();
            }
            solveJunction(prevSpliceSite, exon.getStart() - 1, VariantAnnotationUtils.SPLICE_DONOR_VARIANT,
                    VariantAnnotationUtils.SPLICE_ACCEPTOR_VARIANT, junctionSolution);
            splicing = (splicing || junctionSolution[0]);

            if (variantStart >= exon.getStart()) {
                cdnaExonEnd += (exon.getEnd() - exon.getStart() + 1);
                if (variantStart <= exon.getEnd()) {  // Variant start within the exon
                    cdnaVariantStart = cdnaExonEnd - (exon.getEnd() - variantStart);
                    consequenceType.setCdnaPosition(cdnaVariantStart);
                    variantStartExonNumber = exon.getExonNumber();
                    if (variantEnd <= exon.getEnd()) {  // Both variant start and variant end within the exon  ----||||S|||||E||||----
                        cdnaVariantEnd = cdnaExonEnd - (exon.getEnd() - variantEnd);
                    }
                }
            } else if (variantEnd <= exon.getEnd()) {
                if (variantEnd >= exon.getStart()) {  // Only variant end within the exon  ----||||||||||E||||----
                    cdnaExonEnd += (exon.getEnd() - exon.getStart() + 1);
                    cdnaVariantEnd = cdnaExonEnd - (exon.getEnd() - variantEnd);
                    variantEndExonNumber = exon.getExonNumber();
                } else {  // Variant does not include this exon, variant is located before this exon
                    variantAhead = false;
                }
            } else {  // Variant includes the whole exon. Variant start is located before the exon, variant end is located after the exon
                cdnaExonEnd += (exon.getEnd() - exon.getStart() + 1);
            }
            exonCounter++;
        }
        consequenceType.setExonNumber(variantStartExonNumber != null ? variantStartExonNumber : variantEndExonNumber);
        // Is not intron variant (both ends fall within the same intron)
        if (!junctionSolution[1]) {
            solveExonVariantInPositiveTranscript(splicing, transcriptSequence, cdnaVariantStart, cdnaVariantEnd,
                    firstCdsPhase);
        }

    }

    private void solveNonCodingPositiveTranscript() {
        Exon exon = transcript.getExons().get(0);
        boolean variantAhead = true; // we need a first iteration within the while to ensure junction is solved in case needed
        int cdnaExonEnd = (exon.getEnd() - exon.getStart() + 1);
        int cdnaVariantStart = -1;
        int cdnaVariantEnd = -1;
        int firstCdsPhase = -1;
        boolean[] junctionSolution = {false, false};
        boolean splicing = false;
        Integer variantStartExonNumber = null;
        Integer variantEndExonNumber = null;

        if (transcript.getGenomicCodingStart() <= exon.getEnd()) {
            firstCdsPhase = exon.getPhase();
        }
        if (variantStart >= exon.getStart()) {
            if (variantStart <= exon.getEnd()) {  // Variant start within the exon
                cdnaVariantStart = cdnaExonEnd - (exon.getEnd() - variantStart);
                consequenceType.setCdnaPosition(cdnaVariantStart);
                variantStartExonNumber = exon.getExonNumber();
                if (variantEnd <= exon.getEnd()) {  // Both variant start and variant end within the exon  ----||||S|||||E||||----
                    cdnaVariantEnd = cdnaExonEnd - (exon.getEnd() - variantEnd);
                }
            }
        } else if (variantEnd <= exon.getEnd()) {
            // We do not contemplate that variant end can be located before this exon since this is the first exon
            cdnaVariantEnd = cdnaExonEnd - (exon.getEnd() - variantEnd);
            variantEndExonNumber = exon.getExonNumber();
        } // Variant includes the whole exon. Variant start is located before the exon, variant end is located after the exon


        int exonCounter = 1;
        // This is not a do-while since we cannot call solveJunction until
        while (exonCounter < transcript.getExons().size() && variantAhead) {
            int prevSpliceSite = exon.getEnd() + 1;
            exon = transcript.getExons().get(exonCounter);          // next exon has been loaded
            // Set firsCdsPhase only when the first coding exon is reached
            if (firstCdsPhase == -1 && transcript.getGenomicCodingStart() <= exon.getEnd()) {
                firstCdsPhase = exon.getPhase();
            }
            solveJunction(prevSpliceSite, exon.getStart() - 1, VariantAnnotationUtils.SPLICE_DONOR_VARIANT,
                    VariantAnnotationUtils.SPLICE_ACCEPTOR_VARIANT, junctionSolution);
            splicing = (splicing || junctionSolution[0]);

            if (variantStart >= exon.getStart()) {
                cdnaExonEnd += (exon.getEnd() - exon.getStart() + 1);
                if (variantStart <= exon.getEnd()) {  // Variant start within the exon
                    cdnaVariantStart = cdnaExonEnd - (exon.getEnd() - variantStart);
                    consequenceType.setCdnaPosition(cdnaVariantStart);
                    variantStartExonNumber = exon.getExonNumber();
                    if (variantEnd <= exon.getEnd()) {  // Both variant start and variant end within the exon  ----||||S|||||E||||----
                        cdnaVariantEnd = cdnaExonEnd - (exon.getEnd() - variantEnd);
                    }
                }
            } else if (variantEnd <= exon.getEnd()) {
                if (variantEnd >= exon.getStart()) {  // Only variant end within the exon  ----||||||||||E||||----
                    cdnaExonEnd += (exon.getEnd() - exon.getStart() + 1);
                    cdnaVariantEnd = cdnaExonEnd - (exon.getEnd() - variantEnd);
                    variantEndExonNumber = exon.getExonNumber();
                } else {  // Variant does not include this exon, variant is located before this exon
                    variantAhead = false;
                }
            } else {  // Variant includes the whole exon. Variant start is located before the exon, variant end is located after the exon
                cdnaExonEnd += (exon.getEnd() - exon.getStart() + 1);
            }
            exonCounter++;
        }
        consequenceType.setExonNumber(variantStartExonNumber != null ? variantStartExonNumber : variantEndExonNumber);
        solveMiRNA(cdnaVariantStart, cdnaVariantEnd, junctionSolution[1]);
    }

    private void solveExonVariantInPositiveTranscript(boolean splicing, String transcriptSequence,
                                                      int cdnaVariantStart, int cdnaVariantEnd, int firstCdsPhase) {
        if (variantStart < transcript.getGenomicCodingStart()) {
            // Check transcript has 3 UTR
            if (transcript.getStart() < transcript.getGenomicCodingStart() || transcript.unconfirmedStart()) {
                SoNames.add(VariantAnnotationUtils.FIVE_PRIME_UTR_VARIANT);
            }
            if (variantEnd >= transcript.getGenomicCodingStart()) {
                SoNames.add(VariantAnnotationUtils.CODING_SEQUENCE_VARIANT);
                // cdnaCodingStart<1 if cds_start_NF and phase!=0
                if (transcript.getCdnaCodingStart() > 0 || !transcript.unconfirmedStart()) {
                    SoNames.add(VariantAnnotationUtils.START_LOST);
                }
                if (variantEnd > (transcript.getGenomicCodingEnd() - 3)) {
                    SoNames.add(VariantAnnotationUtils.STOP_LOST);
                    if (variantEnd > transcript.getGenomicCodingEnd()) {
                        // Check transcript has 3 UTR)
                        if (transcript.getEnd() > transcript.getGenomicCodingEnd() || transcript.unconfirmedStart()) {
                            SoNames.add(VariantAnnotationUtils.THREE_PRIME_UTR_VARIANT);
                        }
                    }
                }
            }
        } else if (variantStart <= transcript.getGenomicCodingEnd()) {  // Variant start within coding region
            // Need to define a local cdnaCodingStart because may modified in two lines below
            int cdnaCodingStart = transcript.getCdnaCodingStart();
            if (cdnaVariantStart != -1) {  // cdnaVariantStart may be null if variantStart falls in an intron
                if (transcript.unconfirmedStart()) {
                    cdnaCodingStart -= ((3 - firstCdsPhase) % 3);
                }
                int cdsVariantStart = cdnaVariantStart - cdnaCodingStart + 1;
                consequenceType.setCdsPosition(cdsVariantStart);
                // First place where protein variant annotation is added to the Consequence type,
                // must create the ProteinVariantAnnotation object
                ProteinVariantAnnotation proteinVariantAnnotation = new ProteinVariantAnnotation();
                proteinVariantAnnotation.setPosition(((cdsVariantStart - 1) / 3) + 1);
                consequenceType.setProteinVariantAnnotation(proteinVariantAnnotation);
            }
            if (variantEnd <= transcript.getGenomicCodingEnd()) {  // Variant end also within coding region
                solveCodingExonVariantInPositiveTranscript(splicing, transcriptSequence, cdnaCodingStart,
                        cdnaVariantStart, cdnaVariantEnd);
            } else {
                // Check transcript has 3 UTR)
                if (transcript.getEnd() > transcript.getGenomicCodingEnd() || transcript.unconfirmedEnd()) {
                    SoNames.add(VariantAnnotationUtils.THREE_PRIME_UTR_VARIANT);
                }
                SoNames.add(VariantAnnotationUtils.CODING_SEQUENCE_VARIANT);
                SoNames.add(VariantAnnotationUtils.STOP_LOST);
            }
            // Check transcript has 3 UTR)
        } else if (transcript.getEnd() > transcript.getGenomicCodingEnd() || transcript.unconfirmedEnd()) {
            SoNames.add(VariantAnnotationUtils.THREE_PRIME_UTR_VARIANT);
        }
    }

    private void solveCodingExonVariantInPositiveTranscript(boolean splicing, String transcriptSequence, int cdnaCodingStart,
                                                            int cdnaVariantStart, int cdnaVariantEnd) {
        // This will indicate wether it is needed to add the "coding_sequence_variant" annotation or not
        boolean codingAnnotationAdded = false;

        // cdnaVariantStart=null if variant is intronic. cdnaCodingStart<1 if cds_start_NF and phase!=0
        if (cdnaVariantStart != -1 && cdnaVariantStart < (cdnaCodingStart + 3) && (cdnaCodingStart > 0 || !transcript.unconfirmedStart())) {
            SoNames.add(VariantAnnotationUtils.START_LOST);
            codingAnnotationAdded = true;
        }
        if (cdnaVariantEnd != -1) {
            int finalNtPhase = (transcript.getCdnaCodingEnd() - cdnaCodingStart) % 3;
            Boolean stopToSolve = true;
            // just checks cdnaVariantStart!=null because no splicing means cdnaVariantEnd is also != null
            if (!splicing && cdnaVariantStart != -1) {
                codingAnnotationAdded = true;
                if (Math.abs(variant.getReference().length() - variant.getAlternate().length()) % 3 == 0) {
                    SoNames.add(VariantAnnotationUtils.INFRAME_VARIANT);
                } else {
                    SoNames.add(VariantAnnotationUtils.FRAMESHIFT_VARIANT);
                }
                stopToSolve = false;  // Stop codon annotation will be solved in the line below.
                solveStopCodonPositiveMNV(transcriptSequence, cdnaCodingStart, cdnaVariantStart, cdnaVariantEnd);
            }
            if (cdnaVariantEnd >= (transcript.getCdnaCodingEnd() - finalNtPhase)) {
                if (finalNtPhase != 2) {
//                if (transcript.unconfirmedEnd() && (finalNtPhase != 2)) {
                    SoNames.add(VariantAnnotationUtils.INCOMPLETE_TERMINAL_CODON_VARIANT);
                } else if (stopToSolve) {  // Only if stop codon annotation was not already solved in the if block above
                    SoNames.add(VariantAnnotationUtils.STOP_LOST);
                }
            }
        }
        if (!codingAnnotationAdded) {
            SoNames.add(VariantAnnotationUtils.CODING_SEQUENCE_VARIANT);
        }
    }

    private void solveStopCodonPositiveMNV(String transcriptSequence, int cdnaCodingStart, int cdnaVariantStart,
                                           int cdnaVariantEnd) {
        Integer variantPhaseShift1 = (cdnaVariantStart - cdnaCodingStart) % 3;
        Integer variantPhaseShift2 = (cdnaVariantEnd - cdnaCodingStart) % 3;
        int modifiedCodon1Start = cdnaVariantStart - variantPhaseShift1;
        int modifiedCodon2Start = cdnaVariantEnd - variantPhaseShift2;
        if (modifiedCodon1Start > 0 && (modifiedCodon2Start + 2) <= transcriptSequence.length()) {
            /**
             * Calculate reference codons
             */
            // -1 and +2 because of base 0 String indexing
            String referenceCodon1 = transcriptSequence.substring(modifiedCodon1Start - 1, modifiedCodon1Start + 2);
            String referenceCodon2 = transcriptSequence.substring(modifiedCodon2Start - 1, modifiedCodon2Start + 2);
            int i = cdnaVariantEnd;  // Position (0 based index) in transcriptSequence of the first nt after the deletion
            int codonPosition;

            // Char array to contain the upper/lower-case formatted strings for the codon change, e.g. aGT/ATG
            char[] formattedReferenceCodon1Array = referenceCodon1.toLowerCase().toCharArray();

            /**
             * Calculate codons generated by the insertion
             */
            char[] modifiedCodonArray = referenceCodon1.toCharArray();
            // indexing over transcriptSequence is 0 based, cdnaVariantEnd points to the first position that remains in
            // the allele
            int transcriptSequencePosition = cdnaVariantEnd;
//            int transcriptSequencePosition = cdnaVariantStart;
            int modifiedCodonPosition;
            int modifiedCodonPositionStart = variantPhaseShift1;

            // Char arrays to contain the upper/lower-case formatted strings for the codon change, e.g. aGT/ATG
//            char[] formattedReferenceCodonArray = referenceCodon.toLowerCase().toCharArray();
            char[] formattedModifiedCodonArray = referenceCodon1.toLowerCase().toCharArray();
            boolean useMitochondrialCode = variant.getChromosome().equals("MT");
            boolean firstCodon = true;
            i = 0;

            do {
                for (modifiedCodonPosition = modifiedCodonPositionStart;
                    // Paste alternative nt in the corresponding codon position
                     (modifiedCodonPosition < 3 && i < variant.getAlternate().length()); modifiedCodonPosition++) {
                    modifiedCodonArray[modifiedCodonPosition] = variant.getAlternate().toCharArray()[i];

                    // Edit modified nt to make it upper-case in the formatted strings
                    formattedReferenceCodon1Array[modifiedCodonPosition]
                            = Character.toUpperCase(formattedReferenceCodon1Array[modifiedCodonPosition]);
                    formattedModifiedCodonArray[modifiedCodonPosition]
                            = Character.toUpperCase(variant.getAlternate().toCharArray()[i]);

                    i++;
                }
                for (; modifiedCodonPosition < 3; modifiedCodonPosition++) {  // Concatenate reference codon nts after alternative nts
                    if (transcriptSequencePosition >= transcriptSequence.length()) {
                        int genomicCoordinate = transcript.getEnd() + (transcriptSequencePosition - transcriptSequence.length()) + 1;
                        Query query = new Query(GenomeDBAdaptor.QueryParams.REGION.key(), variant.getChromosome()
                                + ":" + genomicCoordinate
                                + "-" + (genomicCoordinate + 1));
                        modifiedCodonArray[modifiedCodonPosition] = genomeDBAdaptor.getGenomicSequence(query, new QueryOptions())
                                .getResult().get(0).getSequence().charAt(0);
                    } else {
                        modifiedCodonArray[modifiedCodonPosition] = transcriptSequence.charAt(transcriptSequencePosition);
                    }
                    transcriptSequencePosition++;

                    // Edit modified nt to make it upper-case in the formatted strings
                    formattedReferenceCodon1Array[modifiedCodonPosition]
                            = Character.toUpperCase(formattedReferenceCodon1Array[modifiedCodonPosition]);
                    formattedModifiedCodonArray[modifiedCodonPosition]
                            = Character.toUpperCase(modifiedCodonArray[modifiedCodonPosition]);

                }

                // Set codon str, protein ref and protein alt ONLY for the first codon mofified by the insertion
                if (firstCodon) {
                    firstCodon = false;
                    // Only the exact codon where the deletion starts is set
                    consequenceType.setCodon(String.valueOf(formattedReferenceCodon1Array) + "/"
                            + String.valueOf(formattedModifiedCodonArray));
                    // Assumes proteinVariantAnnotation attribute is already initialized
                    consequenceType
                            .getProteinVariantAnnotation()
                            .setReference(VariantAnnotationUtils.getAminoacid(useMitochondrialCode, referenceCodon1));
                    consequenceType
                            .getProteinVariantAnnotation()
                            .setAlternate(VariantAnnotationUtils.getAminoacid(useMitochondrialCode,
                                    String.valueOf(modifiedCodonArray)));
                }

                decideStopCodonModificationAnnotation(SoNames,
                        VariantAnnotationUtils.isStopCodon(useMitochondrialCode, referenceCodon2)
                                ? referenceCodon2 : referenceCodon1, String.valueOf(modifiedCodonArray),
                        useMitochondrialCode);
                modifiedCodonPositionStart = 0;  // Reset the position where the next modified codon must be started to be filled
            } while (i < variant.getAlternate().length());  // All posible new codons generated by the inserted sequence must be checked
        }
    }

    private void solveJunction(Integer spliceSite1, Integer spliceSite2, String leftSpliceSiteTag,
                               String rightSpliceSiteTag, boolean[] junctionSolution) {

        junctionSolution[0] = false;  // Is splicing variant in non-coding region
        junctionSolution[1] = false;  // Variant is intronic and both ends fall within the intron

        // Variant overlaps the rest of intronic region (splice region within the intron and/or rest of intron)
        if (regionsOverlap(spliceSite1 + 2, spliceSite2 - 2, variantStart, variantEnd)) {
            SoNames.add(VariantAnnotationUtils.INTRON_VARIANT);
        }
        if (variantStart >= spliceSite1 && variantEnd <= spliceSite2) {
            junctionSolution[1] = true;  // variant start & end fall within the intron
        }

        if (regionsOverlap(spliceSite1, spliceSite1 + 1, variantStart, variantEnd)) {  // Variant donor/acceptor
            SoNames.add(leftSpliceSiteTag);  // donor/acceptor depending on transcript strand
            // BE CAREFUL: there are introns shorter than 7nts, and even just 1nt long!! (22:36587846)
            junctionSolution[0] = (variantStart <= spliceSite2 || variantEnd <= spliceSite2);
        } else if (regionsOverlap(spliceSite1 + 2, spliceSite1 + 7, variantStart, variantEnd)) {
            SoNames.add(VariantAnnotationUtils.SPLICE_REGION_VARIANT);
        } else if (regionsOverlap(spliceSite1 - 3, spliceSite1 - 1, variantStart, variantEnd)) {
            SoNames.add(VariantAnnotationUtils.SPLICE_REGION_VARIANT);
        }

        if (regionsOverlap(spliceSite2 - 1, spliceSite2, variantStart, variantEnd)) {  // Variant donor/acceptor
            SoNames.add(rightSpliceSiteTag);  // donor/acceptor depending on transcript strand
            // BE CAREFUL: there are introns shorter than 7nts, and even just 1nt long!! (22:36587846)
            junctionSolution[0] = (spliceSite1 <= variantStart || spliceSite1 <= variantEnd);
        } else if (regionsOverlap(spliceSite2 - 7, spliceSite2 - 2, variantStart, variantEnd)) {
            SoNames.add(VariantAnnotationUtils.SPLICE_REGION_VARIANT);
            // BE CAREFUL: there are introns shorter than 7nts, and even just 1nt long!! (22:36587846)
            junctionSolution[0] = (spliceSite1 <= variantStart || spliceSite1 <= variantEnd);
        } else if (regionsOverlap(spliceSite2 + 1, spliceSite2 + 3, variantStart, variantEnd)) {
            SoNames.add(VariantAnnotationUtils.SPLICE_REGION_VARIANT);
        }
    }

}
