# GO-enrichment hero candidates — per-clade shortlist

Source: `out/significant_terms_annotated.tsv` (post-rebuild, 28 clades).
Gate: `fold >= 3`, namespace in {BP, MF, CC}, sigma=2 on the underlying
`unique_pairs` filter. Top 3 per clade ranked by `hero_score = −log10(q) × log2(fold) × √k`,
preferring clade-specific terms (term only reaches q ≤ 0.01 in one clade).

## Top-5 narrative picks

1. **Bilateria — mitochondrial Ca²⁺ import** (`GO:0036444`, fold 71, q 1.4×10⁻³). The mitochondrial calcium uniporter (MCU) complex is textbook bilaterian-specific — our method recovers it cold. Strong proof-of-approach hit.
2. **Ctenophora — transporter complex** (`GO:1990351`, fold 163, q 8×10⁻⁴) and **H3K36 demethylase** (`GO:0051864`, fold 76, q 2×10⁻³). Largest fold enrichment in the panel; fits ctenophore-unique neurogenesis/physiology literature.
3. **Platyhelminthes — apoptotic regulation + synapse pruning** (`GO:0042981`, k=4, q 1×10⁻⁵, fold 26; `GO:1905806`, fold 120, q 1×10⁻⁴). Flatworm regenerative/neural-plasticity narrative.
4. **Cnidaria — angiogenesis/cell-migration regulation** (`GO:0045765`, fold 53; `GO:0030334`, fold 30). Deep-conserved morphogenetic machinery in a clade without blood vessels — evo-devo story.
5. **Nematoda — minus-end microtubule motor + nucleotide binding** (`GO:0008569`, fold 90, q 4×10⁻⁵; `GO:0000166`, k=8, fold 12, q 5×10⁻⁶). Strongest q in the panel; ties to nematode spermatogenesis / sensory cilia.

## Full per-clade shortlist (top 3 per clade)

| clade | rank | go_id | ns | k | fold | q | clade-specific | go_name |
|---|---|---|---|---|---|---|---|---|
| Annelida | 1 | `GO:0033186` | CC | 2 | 80.3 | 2.92e-03 | ✓ | CAF-1 complex |
| Annelida | 2 | `GO:0015078` | MF | 2 | 34.4 | 9.36e-03 | ✓ | proton transmembrane transporter activity |
| Annelida | 3 | `GO:0071339` | CC | 2 | 15.0 | 5.10e-02 | ✓ | MLL1 complex |
| Arthropoda | 1 | `GO:0062072` | MF | 2 | 24.1 | 9.36e-03 | ✓ | histone H3K9me2/3 reader activity |
| Arthropoda | 2 | `GO:0008047` | MF | 3 | 12.1 | 9.36e-03 | ✓ | enzyme activator activity |
| Arthropoda | 3 | `GO:0045505` | MF | 6 | 5.8 | 1.74e-02 | ✓ | dynein intermediate chain binding |
| Bilateria | 1 | `GO:0036444` | BP | 2 | 71.6 | 1.36e-03 | ✓ | calcium import into the mitochondrion |
| Bilateria | 2 | `GO:0019855` | MF | 2 | 71.6 | 2.72e-03 | ✓ | calcium channel inhibitor activity |
| Bilateria | 3 | `GO:0003281` | BP | 2 | 50.9 | 6.42e-03 | ✓ | ventricular septum development |
| Bivalvia | 1 | `GO:0016251` | MF | 2 | 47.7 | 7.62e-03 | ✓ | RNA polymerase II general transcription initiation factor activity |
| Bivalvia | 2 | `GO:0000151` | CC | 3 | 18.5 | 5.92e-03 | ✓ | ubiquitin ligase complex |
| Bivalvia | 3 | `GO:0005669` | CC | 2 | 34.7 | 8.24e-03 | ✓ | transcription factor TFIID complex |
| Cephalopoda | 1 | `GO:0045190` | BP | 2 | 21.8 | 2.06e-02 | ✓ | isotype switching |
| Cephalopoda | 2 | `GO:0006979` | BP | 2 | 11.3 | 1.28e-02 | ✓ | response to oxidative stress |
| Cephalopoda | 3 | `GO:0002181` | BP | 5 | 4.4 | 3.14e-02 | ✓ | cytoplasmic translation |
| Chordata | 1 | `GO:0060173` | BP | 2 | 29.9 | 3.02e-03 | ✓ | limb development |
| Chordata | 2 | `GO:0034063` | BP | 2 | 32.7 | 5.19e-03 | ✓ | stress granule assembly |
| Chordata | 3 | `GO:0006325` | BP | 3 | 12.3 | 2.80e-03 | ✓ | chromatin organization |
| Clitellata | 1 | `GO:0019477` | BP | 2 | 60.3 | 1.55e-03 | ✓ | L-lysine catabolic process |
| Clitellata | 2 | `GO:0032956` | BP | 2 | 22.0 | 7.78e-03 | ✓ | regulation of actin cytoskeleton organization |
| Clitellata | 3 | `GO:0060173` | BP | 2 | 11.6 | 1.25e-01 | ✓ | limb development |
| Cnidaria | 1 | `GO:0045765` | BP | 2 | 53.3 | 7.58e-03 | ✓ | regulation of angiogenesis |
| Cnidaria | 2 | `GO:0030334` | BP | 2 | 30.5 | 3.50e-03 | ✓ | regulation of cell migration |
| Cnidaria | 3 | `GO:1904781` | BP | 2 | 50.9 | 7.55e-03 | ✓ | positive regulation of protein localization to centrosome |
| Coleoidea | 1 | `GO:0006979` | BP | 2 | 15.8 | 6.55e-03 | ✓ | response to oxidative stress |
| Coleoidea | 2 | `GO:0045190` | BP | 2 | 13.7 | 4.08e-02 | ✓ | isotype switching |
| Coleoidea | 3 | `GO:0002181` | BP | 4 | 4.6 | 4.71e-02 | ✓ | cytoplasmic translation |
| Ctenophora | 1 | `GO:1990351` | CC | 2 | 163.6 | 7.99e-04 | ✓ | transporter complex |
| Ctenophora | 2 | `GO:0051864` | MF | 2 | 76.3 | 1.95e-03 | ✓ | histone H3K36 demethylase activity |
| Ctenophora | 3 | `GO:1990531` | CC | 2 | 76.3 | 2.47e-03 | ✓ | phospholipid-translocating ATPase complex |
| Decapodiformes | 1 | `GO:0071817` | CC | 2 | 114.5 | 1.23e-03 | ✓ | MMXD complex |
| Decapodiformes | 2 | `GO:0007059` | BP | 2 | 16.4 | 9.29e-03 | ✓ | chromosome segregation |
| Decapodiformes | 3 | `GO:0006979` | BP | 4 | 5.6 | 2.83e-02 | ✓ | response to oxidative stress |
| Deuterostomia | 1 | `GO:0005852` | CC | 6 | 8.9 | 7.90e-04 | ✓ | eukaryotic translation initiation factor 3 complex |
| Deuterostomia | 2 | `GO:0033290` | CC | 6 | 8.9 | 7.90e-04 | ✓ | eukaryotic 48S preinitiation complex |
| Deuterostomia | 3 | `GO:0001568` | BP | 2 | 58.7 | 1.30e-03 | ✓ | blood vessel development |
| Diptera | 1 | `GO:0033762` | BP | 2 | 44.9 | 1.36e-02 | ✓ | response to glucagon |
| Diptera | 2 | `GO:0007059` | BP | 2 | 18.2 | 1.26e-02 | ✓ | chromosome segregation |
| Diptera | 3 | `GO:0051864` | MF | 2 | 29.9 | 4.17e-02 | ✓ | histone H3K36 demethylase activity |
| Echinodermata | 1 | `GO:0003755` | MF | 2 | 83.3 | 6.25e-04 | ✓ | peptidyl-prolyl cis-trans isomerase activity |
| Echinodermata | 2 | `GO:0031146` | BP | 2 | 22.0 | 1.37e-02 | ✓ | SCF-dependent proteasomal ubiquitin-dependent protein catabolic process |
| Echinodermata | 3 | `GO:0002181` | BP | 2 | 9.0 | 4.79e-02 | ✓ | cytoplasmic translation |
| Gastropoda | 1 | `GO:0031167` | BP | 2 | 81.8 | 4.14e-04 | ✓ | rRNA methylation |
| Gastropoda | 2 | `GO:0006851` | BP | 4 | 17.2 | 1.19e-03 | ✓ | mitochondrial calcium ion transmembrane transport |
| Gastropoda | 3 | `GO:0006814` | BP | 2 | 45.8 | 2.14e-03 | ✓ | sodium ion transport |
| Hexapoda | 1 | `GO:0072711` | BP | 2 | 52.6 | 5.30e-03 | ✓ | cellular response to hydroxyurea |
| Hexapoda | 2 | `GO:0005802` | CC | 3 | 8.7 | 1.41e-02 | ✓ | trans-Golgi network |
| Hexapoda | 3 | `GO:0008289` | MF | 2 | 14.1 | 1.49e-02 | ✓ | lipid binding |
| Insecta | 1 | `GO:0042147` | BP | 2 | 13.6 | 8.87e-03 | ✓ | retrograde transport, endosome to Golgi |
| Insecta | 2 | `GO:0003777` | MF | 9 | 3.6 | 1.48e-02 | ✓ | microtubule motor activity |
| Insecta | 3 | `GO:0005802` | CC | 3 | 8.7 | 1.41e-02 | ✓ | trans-Golgi network |
| Mammalia | 1 | `GO:0006515` | BP | 2 | 57.2 | 4.30e-04 | ✓ | protein quality control for misfolded or incompletely synthesized proteins |
| Mammalia | 2 | `GO:0004553` | MF | 2 | 48.2 | 5.14e-03 | ✓ | hydrolase activity, hydrolyzing O-glycosyl compounds |
| Mammalia | 3 | `GO:0055088` | BP | 2 | 31.2 | 4.09e-03 | ✓ | lipid homeostasis |
| Mollusca | 1 | `GO:0071339` | CC | 2 | 58.7 | 6.61e-03 | ✓ | MLL1 complex |
| Mollusca | 2 | `GO:0045893` | BP | 6 | 6.2 | 3.19e-03 | ✓ | positive regulation of DNA-templated transcription |
| Mollusca | 3 | `GO:0034446` | BP | 2 | 30.1 | 7.52e-03 | ✓ | substrate adhesion-dependent cell spreading |
| Nematoda | 1 | `GO:0000166` | MF | 8 | 12.3 | 4.59e-06 | ✓ | nucleotide binding |
| Nematoda | 2 | `GO:0008569` | MF | 3 | 90.4 | 4.24e-05 | ✓ | minus-end-directed microtubule motor activity |
| Nematoda | 3 | `GO:0030286` | CC | 3 | 72.3 | 1.15e-04 | ✓ | dynein complex |
| Neoptera | 1 | `GO:0045780` | BP | 2 | 61.9 | 3.56e-03 | ✓ | positive regulation of bone resorption |
| Neoptera | 2 | `GO:0090090` | BP | 2 | 17.0 | 8.65e-03 | ✓ | negative regulation of canonical Wnt signaling pathway |
| Neoptera | 3 | `GO:0042147` | BP | 2 | 12.7 | 1.02e-02 | ✓ | retrograde transport, endosome to Golgi |
| Platyhelminthes | 1 | `GO:0042981` | BP | 4 | 26.3 | 9.59e-06 | ✓ | regulation of apoptotic process |
| Platyhelminthes | 2 | `GO:1905806` | BP | 2 | 120.5 | 1.06e-04 | ✓ | regulation of synapse pruning |
| Platyhelminthes | 3 | `GO:0070585` | BP | 2 | 26.8 | 9.07e-03 | ✓ | protein localization to mitochondrion |
| Porifera | 1 | `GO:0070778` | BP | 2 | 120.5 | 2.22e-04 | ✓ | L-aspartate transmembrane transport |
| Porifera | 2 | `GO:0015810` | BP | 2 | 120.5 | 2.22e-04 | ✓ | aspartate transmembrane transport |
| Porifera | 3 | `GO:0015813` | BP | 2 | 120.5 | 2.22e-04 | ✓ | L-glutamate transmembrane transport |
| Protostomia | 1 | `GO:2000138` | BP | 2 | 61.9 | 1.84e-03 | ✓ | positive regulation of cell proliferation involved in heart morphogenesis |
| Protostomia | 2 | `GO:0060716` | BP | 2 | 61.9 | 1.84e-03 | ✓ | labyrinthine layer blood vessel development |
| Protostomia | 3 | `GO:0036302` | BP | 2 | 41.3 | 3.13e-03 | ✓ | atrioventricular canal development |
| Scaphopoda | 1 | `GO:0001535` | CC | 2 | 42.4 | 9.01e-03 | ✓ | radial spoke head |
| Scaphopoda | 2 | `GO:0160112` | CC | 5 | 9.2 | 6.90e-03 | ✓ | axonemal B tubule inner sheath |
| Scaphopoda | 3 | `GO:0006914` | BP | 3 | 17.0 | 8.48e-03 | ✓ | autophagy |
| Spiralia | 1 | `GO:0065003` | BP | 2 | 23.4 | 8.75e-03 | ✓ | protein-containing complex assembly |
| Spiralia | 2 | `GO:0000151` | CC | 2 | 21.1 | 1.52e-02 | ✓ | ubiquitin ligase complex |
| Spiralia | 3 | `GO:2000481` | BP | 2 | 11.0 | 1.11e-01 | ✓ | positive regulation of cAMP-dependent protein kinase activity |
| Teleostei | 1 | `GO:0008757` | MF | 3 | 68.7 | 7.29e-05 | ✓ | S-adenosylmethionine-dependent methyltransferase activity |
| Teleostei | 2 | `GO:0008173` | MF | 3 | 57.2 | 7.29e-05 | ✓ | RNA methyltransferase activity |
| Teleostei | 3 | `GO:0008168` | MF | 3 | 18.6 | 2.91e-03 | ✓ | methyltransferase activity |
| Vertebrata | 1 | `GO:0008344` | BP | 2 | 42.4 | 1.99e-03 | ✓ | adult locomotory behavior |
| Vertebrata | 2 | `GO:0010977` | BP | 2 | 31.8 | 1.99e-03 | ✓ | negative regulation of neuron projection development |
| Vertebrata | 3 | `GO:0032456` | BP | 8 | 5.8 | 5.22e-03 | ✓ | endocytic recycling |
