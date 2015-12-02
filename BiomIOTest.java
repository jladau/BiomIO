package edu.ucsf.io;

import static org.junit.Assert.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import org.junit.Test;

import edu.ucsf.io.BiomIO;



public class BiomIOTest {

	/**Directory with test data sets.**/
	private String sTestDataDir;
	
	/**Test data file.**/
	private String sTestFile;
	
	/**Full table. Assumed not to be rarefied, with OTUs at multiple taxonomic resolutions, **/
	private BiomIO bio1;
	
	/**rgsObservations = correct set of observation IDs after filtering**/
	private String rgsCorrectObservationIDs[];
	
	/**rgsSamples = correct set of sample IDs after filtering**/
	private String rgsCorrectSampleIDs[];
	
	/**rgdData = correct set of data after filtering**/
	private double rgdCorrectData[][];
	
	/**
	 * Constructor.
	 * @param sTestDataDir Directory with test data sets.
	 */
	public BiomIOTest(){
		sTestDataDir="/home/jladau/Documents/Research/Data/Microbial_Community_Samples/BiomIOTestData";
		sTestFile="/rich_sparse_otu_table_hdf5.biom";
		bio1=new BiomIO(sTestDataDir + "/" + sTestFile);
	}

	@Test
	public void checkRarefied_TableIsRarefied_ReturnsTrue(){
		try {
			bio1.rarefy(3);
		} catch (Exception e) {
			fail(e.getMessage());
		}
		assertTrue(bio1.checkRarefied());
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
	}
	
	@Test
	public void checkRarefied_TableIsNotRarefied_ReturnsFalse(){
		assertFalse(bio1.checkRarefied());
	}
	
	@Test
	public void resampleWithReplacement_TableIsResampled_TableIsCorrect(){
		
		bio1 = bio1.resampleWithReplacement(1234);
		rgsCorrectSampleIDs=new String[]{"Sample3.1","Sample6.1","Sample6.2","Sample1.1","Sample3.2","Sample4.1"};
		rgsCorrectObservationIDs=new String[]{"GG_OTU_1","GG_OTU_2","GG_OTU_3","GG_OTU_4","GG_OTU_5"};
		rgdCorrectData = new double[][]{
				{1,0,0,0,1,0},
				{0,1,1,5,0,2},
				{1,2,2,0,1,4},
				{1,1,1,2,1,0},
				{1,0,0,0,1,0}};
		this.checkTableIsCorrect();
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
	}
	
	@Test
	public void collapse_TableIsCollapsedByTaxon_TableIsCorrect(){
		
		bio1.collapse("kingdom", bio1.axsObservation, false);
		rgsCorrectObservationIDs = new String[]{"k__Bacteria", "k__Archaea"};
		rgsCorrectSampleIDs=new String[]{"Sample1","Sample2","Sample3","Sample4","Sample5","Sample6"};
		rgdCorrectData=new double[][]{{7,3,3,2,3,2},{0,0,1,4,0,2}};
		this.checkTableIsCorrect();
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
	}
	
	@Test
	public void convertToPresenceAbsence_TableIsPresenceAbsence_TableIsCorrect(){
		
		bio1.convertToPresenceAbsence();
		rgsCorrectSampleIDs=new String[]{"Sample1","Sample2","Sample3","Sample4","Sample5","Sample6"};
		rgsCorrectObservationIDs=new String[]{"GG_OTU_1","GG_OTU_2","GG_OTU_3","GG_OTU_4","GG_OTU_5"};
		rgdCorrectData = new double[][]{
				{0,0,1,0,0,0},
				{1,1,0,1,1,1},
				{0,0,1,1,0,1},
				{1,1,1,0,0,1},
				{0,1,1,0,0,0}};		
		this.checkTableIsCorrect();
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
	}
	
	@Test
	public void getRichness_RichnessGotten_RichnessCorrect(){
		
		//map1 = map from function
		//map2 = correct map
		
		HashMap<String,Double> map1;
		HashMap<String,Double> map2;
		
		map1 = new HashMap<String,Double>();
		map1.put("Sample1", 2.);
		map1.put("Sample2", 3.);
		map1.put("Sample3", 4.);
		map1.put("Sample4", 2.);
		map1.put("Sample5", 1.);
		map1.put("Sample6", 3.);
		
		map2 = bio1.getRichness();
		
		for(String s:map2.keySet()){
			assertEquals(map1.get(s),map2.get(s),0.0000000001);
		}
	}
	
	@Test
	public void equals_TablesAreEquals_ReturnsTrue(){
		assertTrue(bio1.equals(bio1));
	}
	
	@Test
	public void filter_SamplesAreFiltered_TableIsCorrect(){
		
		//set1 = set of samples to filter by
		
		HashSet<String> set1;
		
		rgsCorrectSampleIDs=new String[]{"Sample1","Sample3","Sample5"};
		rgsCorrectObservationIDs=new String[]{"GG_OTU_1","GG_OTU_2","GG_OTU_3","GG_OTU_4","GG_OTU_5"};
		rgdCorrectData = new double[][]{
				{0,1,0},
				{5,0,3},
				{0,1,0},
				{2,1,0},
				{0,1,0}};
		set1 = new HashSet<String>();
		for(int i=0;i<rgsCorrectSampleIDs.length;i++){
			set1.add(rgsCorrectSampleIDs[i]);
		}
		try{
			bio1.filter(set1, bio1.axsSample);
		}catch(Exception e){
			fail(e.getMessage());
		}
		this.checkTableIsCorrect();
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
	}
	
	@Test
	public void filter_ObservationsAreFiltered_TableIsCorrect(){
		
		//set1 = set of observations to filter by
		
		HashSet<String> set1;
		
		rgsCorrectSampleIDs=new String[]{"Sample1","Sample2","Sample3","Sample4","Sample5","Sample6"};
		rgsCorrectObservationIDs=new String[]{"GG_OTU_1","GG_OTU_3","GG_OTU_4"};
		rgdCorrectData = new double[][]{
				{0,0,1,0,0,0},
				{0,0,1,4,0,2},
				{2,1,1,0,0,1}};
		set1 = new HashSet<String>();
		for(int i=0;i<rgsCorrectObservationIDs.length;i++){
			set1.add(rgsCorrectObservationIDs[i]);
		}
		try{
			bio1.filter(set1, bio1.axsObservation);
		}catch(Exception e){
			fail(e.getMessage());
		}
		this.checkTableIsCorrect();
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
	}
	
	@Test
	public void filterByNoMetadata_SamplesAreFiltered_TableIsCorrect(){
		
		rgsCorrectSampleIDs=new String[]{"Sample1","Sample2","Sample3","Sample6"};
		rgsCorrectObservationIDs=new String[]{"GG_OTU_1","GG_OTU_2","GG_OTU_3","GG_OTU_4","GG_OTU_5"};
		rgdCorrectData = new double[][]{
				{0,0,1,0},
				{5,1,0,1},
				{0,0,1,2},
				{2,1,1,1},
				{0,1,1,0}};
		try{
			bio1.filterByNoMetadata(new String[]{"LinkerPrimerSequence","Description","BarcodeSequence"}, bio1.axsSample);
		}catch(Exception e){
			fail(e.getMessage());
		}
		this.checkTableIsCorrect();
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
	}
	
	@Test
	public void filterByPrevalence_ObservationsAreFiltered_TableIsCorrect(){
		
		rgsCorrectSampleIDs=new String[]{"Sample1","Sample2","Sample3","Sample4","Sample5","Sample6"};
		rgsCorrectObservationIDs=new String[]{"GG_OTU_2","GG_OTU_3","GG_OTU_4"};
		rgdCorrectData = new double[][]{
				{5,1,0,2,3,1},
				{0,0,1,4,0,2},
				{2,1,1,0,0,1}};
		try{
			bio1.filterByPrevalence(3);
		}catch(Exception e){
			fail(e.getMessage());
		}
		this.checkTableIsCorrect();
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
	}
	
	@Test
	public void filterFromFile_SamplesAreFiltered_SampleIDsAreCorrect(){
		
		//lstOut = current filter output
		
		ArrayList<String> lstOut;
		
		rgsCorrectSampleIDs=new String[]{"Sample1","Sample2","Sample3"};
		rgsCorrectObservationIDs=new String[]{"GG_OTU_1","GG_OTU_2","GG_OTU_3","GG_OTU_4","GG_OTU_5"};
		rgdCorrectData = new double[][]{
				{0,0,1},
				{5,1,0},
				{0,0,1},
				{2,1,1},
				{0,1,1}};
		try{
			lstOut = new ArrayList<String>();
			lstOut.add("Sample1");
			lstOut.add("Sample2");
			lstOut.add("Sample3");
			DataIO.writeToFile(lstOut, "/tmp/SampleFilter.txt");
			bio1.filterFromFile("/tmp/SampleFilter.txt", bio1.axsSample);
		}catch(Exception e){
			fail(e.getMessage());
		}
		this.checkTableIsCorrect();
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
	}
	
	@Test
	public void filterFromFile_ObservationsAreFiltered_MatrixIsCorrect(){

		//lstOut = current filter output
		
		ArrayList<String> lstOut;
		
		rgsCorrectSampleIDs=new String[]{"Sample1","Sample2","Sample3","Sample4","Sample5","Sample6"};
		rgsCorrectObservationIDs=new String[]{"GG_OTU_2","GG_OTU_3","GG_OTU_4"};
		rgdCorrectData = new double[][]{
				{5,1,0,2,3,1},
				{0,0,1,4,0,2},
				{2,1,1,0,0,1}};
		try{
			lstOut = new ArrayList<String>();
			lstOut.add("GG_OTU_2");
			lstOut.add("GG_OTU_3");
			lstOut.add("GG_OTU_4");
			DataIO.writeToFile(lstOut, "/tmp/ObservationFilter.txt");
			bio1.filterFromFile("/tmp/ObservationFilter.txt", bio1.axsObservation);
		}catch(Exception e){
			fail(e.getMessage());
		}
		this.checkTableIsCorrect();
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
	}

	@Test
	public void getItem_ObservationIsGotten_RowIsCorrect(){
		
		//map1 = correct row
		
		HashMap<String,Double> map1;
		
		map1 = new HashMap<String,Double>();
		map1.put("Sample1", 5.);
		map1.put("Sample2", 1.);
		map1.put("Sample3", 0.);
		map1.put("Sample4", 2.);
		map1.put("Sample5", 3.);
		map1.put("Sample6", 1.);
		this.checkTableMapIsCorrect(map1, bio1.getItem(bio1.axsObservation, "GG_OTU_2"));
	}
	
	@Test
	public void getItem_SampleIsGotten_ColumnIsCorrect(){
		
		//map1 = correct row
		
		HashMap<String,Double> map1;
		
		map1 = new HashMap<String,Double>();
		map1.put("GG_OTU_1", 0.);
		map1.put("GG_OTU_2", 2.);
		map1.put("GG_OTU_3", 4.);
		map1.put("GG_OTU_4", 0.);
		map1.put("GG_OTU_5", 0.);
		this.checkTableMapIsCorrect(map1, bio1.getItem(bio1.axsSample, "Sample4"));
	}
	
	@Test
	public void getMean_ObservationMeanIsTaken_ValueIsCorrect(){
		assertEquals(2.,bio1.getMean(bio1.axsObservation, "GG_OTU_2"),0.00000001);
	}
	
	@Test
	public void getMean_SampleMeanIsTaken_ValueIsCorrect(){
		assertEquals(1.2,bio1.getMean(bio1.axsSample, "Sample4"),0.00000001);
	}
	
	@Test
	public void getMeans_ObservationMeansAreTaken_ValuesAreCorrect(){
		
		//map1 = correct means
		
		HashMap<String,Double> map1;
		
		map1 = new HashMap<String,Double>();
		map1.put("GG_OTU_1", 0.16666666666667);
		map1.put("GG_OTU_2", 2.);
		map1.put("GG_OTU_3", 1.16666666666667);
		map1.put("GG_OTU_4", 0.83333333333333);
		map1.put("GG_OTU_5", 0.33333333333333);
		this.checkTableMapIsCorrect(map1, bio1.getMeans(bio1.axsObservation));
	}
	
	@Test
	public void getMeans_SampleMeansAreTaken_ValuesAreCorrect(){
		
		//map1 = correct means
		
		HashMap<String,Double> map1;
		
		map1 = new HashMap<String,Double>();
		map1.put("Sample1", 1.4);
		map1.put("Sample2", 0.6);
		map1.put("Sample3", 0.8);
		map1.put("Sample4", 1.2);
		map1.put("Sample5", 0.6);
		map1.put("Sample6", 0.8);
		
		this.checkTableMapIsCorrect(map1, bio1.getMeans(bio1.axsSample));
	}
	
	@Test
	public void getNonzeroCount_ObservationCountIsTaken_ValueIsCorrect(){
		assertEquals(5.,bio1.getNonzeroCount(bio1.axsObservation, "GG_OTU_2"),0.00000001);
	}
	
	@Test
	public void getNonzeroCount_SampleCountIsTaken_ValueIsCorrect(){
		assertEquals(2,bio1.getNonzeroCount(bio1.axsSample, "Sample4"),0.00000001);
	}
	
	@Test
	public void getNonzeroCounts_ObservationCountsAreTaken_ValuesAreCorrect(){
		
		//map1 = correct counts
		//map2 = observed counts (double)
		//map3 = observted counts (int)
		
		HashMap<String,Integer> map3;
		HashMap<String,Double> map2;
		HashMap<String,Double> map1;
		
		map1 = new HashMap<String,Double>();
		map1.put("GG_OTU_1", 1.);
		map1.put("GG_OTU_2", 5.);
		map1.put("GG_OTU_3", 3.);
		map1.put("GG_OTU_4", 4.);
		map1.put("GG_OTU_5", 2.);
		map2 = new HashMap<String,Double>();
		map3 = bio1.getNonzeroCounts(bio1.axsObservation);
		for(String s:map3.keySet()){
			map2.put(s, (double) map3.get(s));
		}
		this.checkTableMapIsCorrect(map1, map2);
	}
	
	@Test
	public void getNonzeroCounts_SampleCountsAreTaken_ValuesAreCorrect(){
		
		//map1 = correct counts
		//map2 = observed counts (double)
		//map3 = observted counts (int)
		
		HashMap<String,Integer> map3;
		HashMap<String,Double> map2;
		HashMap<String,Double> map1;
		
		map1 = new HashMap<String,Double>();
		map1.put("Sample1", 2.);
		map1.put("Sample2", 3.);
		map1.put("Sample3", 4.);
		map1.put("Sample4", 2.);
		map1.put("Sample5", 1.);
		map1.put("Sample6", 3.);
		map2 = new HashMap<String,Double>();
		map3 = bio1.getNonzeroCounts(bio1.axsSample);
		for(String s:map3.keySet()){
			map2.put(s, (double) map3.get(s));
		}
		this.checkTableMapIsCorrect(map1, map2);
	}
	
	@Test
	public void getShape_IsRun_ValueIsCorrect(){
		assertArrayEquals(new int[]{5,6},bio1.getShape());
	}
	
	@Test
	public void getValueByIDs_RandomValuesSelected_ValuesAreCorrect(){
		assertEquals(2.,bio1.getValueByIDs("GG_OTU_2", "Sample4"),0.00001);
		assertEquals(4.,bio1.getValueByIDs("GG_OTU_3", "Sample4"),0.00001);
		assertEquals(0.,bio1.getValueByIDs("GG_OTU_5", "Sample1"),0.00001);
		assertEquals(1.,bio1.getValueByIDs("GG_OTU_5", "Sample3"),0.00001);
		assertEquals(0.,bio1.getValueByIDs("GG_OTU_1", "Sample6"),0.00001);
	}
	
	@Test
	public void getValueByIndices_RandomValuesSelected_ValuesAreCorrect(){
		assertEquals(2.,bio1.getValueByIndices(1,3),0.00001);
		assertEquals(4.,bio1.getValueByIndices(2,3),0.00001);
		assertEquals(0.,bio1.getValueByIndices(4,0),0.00001);
		assertEquals(1.,bio1.getValueByIndices(4,1),0.00001);
		assertEquals(0.,bio1.getValueByIndices(0,5),0.00001);
	}
	
	@Test
	public void normalize_Normalized_TableIsCorrect(){
		rgsCorrectSampleIDs=new String[]{"Sample1","Sample2","Sample3","Sample4","Sample5","Sample6"};
		rgsCorrectObservationIDs=new String[]{"GG_OTU_1","GG_OTU_2","GG_OTU_3","GG_OTU_4","GG_OTU_5"};
		rgdCorrectData = new double[][]{
				{0,0,0.25,0,0,0},
				{0.714285714285714,0.333333333333333,0,0.333333333333333,1,0.25},
				{0,0,0.25,0.666666666666667,0,0.5},
				{0.285714285714286,0.333333333333333,0.25,0,0,0.25},
				{0,0.333333333333333,0.25,0,0,0}};
		bio1.normalize();
		this.checkTableIsCorrect();
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
	}
	
	@Test
	public void getShannon_ShannonGotten_ShannonCorrect(){
		
		//map1 = map from function
		//map2 = correct map
		
		HashMap<String,Double> map1;
		HashMap<String,Double> map2;
		
		map1 = new HashMap<String,Double>();
		map1.put("Sample1", 0.59826958858526);
		map1.put("Sample2", 1.09861228866811);
		map1.put("Sample3", 1.38629436111989);
		map1.put("Sample4", 0.63651416829481);
		map1.put("Sample5", 0.);
		map1.put("Sample6", 1.03972077083992);
		
		bio1.normalize();
		map2 = bio1.getShannon();
		for(String s:map2.keySet()){
			assertEquals(map1.get(s),map2.get(s),0.0000000001);
		}
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
	}
	
	//TODO add test for printMetadata
	
	@Test
	public void printTable_TableIsPrinted_OutputIsCorrect(){
		
		//lst1 = correct output
		//lst2 = output
		
		ArrayList<String> lst1;
		ArrayList<String> lst2;
		
		lst1 = new ArrayList<String>();
		lst1.add("# Constructed from biom file");
		lst1.add("#OTU ID,Sample1,Sample2,Sample3,Sample4,Sample5,Sample6");
		lst1.add("GG_OTU_1,0.0,0.0,1.0,0.0,0.0,0.0");
		lst1.add("GG_OTU_2,5.0,1.0,0.0,2.0,3.0,1.0");
		lst1.add("GG_OTU_3,0.0,0.0,1.0,4.0,0.0,2.0");
		lst1.add("GG_OTU_4,2.0,1.0,1.0,0.0,0.0,1.0");
		lst1.add("GG_OTU_5,0.0,1.0,1.0,0.0,0.0,0.0");
		lst2 = bio1.printTable();
		for(int i=0;i<lst1.size();i++){
			assertEquals(lst1.get(i),lst2.get(i));
		}
	}
	
	@Test
	public void rarefy_TableIsRarefied_TableIsCorrect(){
		
		//mapSum = sums
		
		HashMap<String,Double> mapSum;
		
		try {
			bio1.rarefy(5);
		} catch (Exception e) {
			fail(e.getMessage());
		}
		rgsCorrectSampleIDs=new String[]{"Sample1","Sample4"};
		this.checkAxisIDsAreCorrect(rgsCorrectSampleIDs, bio1.axsSample);
		mapSum = bio1.sum(bio1.axsSample);
		for(String s:mapSum.keySet()){
			assertEquals(5.,mapSum.get(s),0.0000000001);
		}
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
	}
	
	@Test
	public void subsample_TableIsSubsampled_TableIsCorrect(){
		
		//mapSum = sums
		
		HashMap<String,Double> mapSum;
		
		try {
			bio1.subsample(6, bio1.axsObservation);
		} catch (Exception e) {
			fail(e.getMessage());
		}
		rgsCorrectObservationIDs=new String[]{"GG_OTU_2","GG_OTU_3"};
		this.checkAxisIDsAreCorrect(rgsCorrectObservationIDs, bio1.axsObservation);
		mapSum = bio1.sum(bio1.axsObservation);
		for(String s:mapSum.keySet()){
			assertEquals(6.,mapSum.get(s),0.0000000001);
		}
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
	}
	
	@Test
	public void sum_TableIsSummed_SumsAreCorrect(){
		
		//mapSum = sums
		//mapSumCorrect = correct sums
		
		HashMap<String,Double> mapSum;
		HashMap<String,Double> mapSumCorrect;
		
		mapSumCorrect = new HashMap<String,Double>();
		mapSumCorrect.put("GG_OTU_1", 1.);
		mapSumCorrect.put("GG_OTU_2", 12.);
		mapSumCorrect.put("GG_OTU_3", 7.);
		mapSumCorrect.put("GG_OTU_4", 5.);
		mapSumCorrect.put("GG_OTU_5", 2.);
		mapSumCorrect.put("Sample1", 7.);
		mapSumCorrect.put("Sample2", 3.);
		mapSumCorrect.put("Sample3", 4.);
		mapSumCorrect.put("Sample4", 6.);
		mapSumCorrect.put("Sample5", 3.);
		mapSumCorrect.put("Sample6", 4.);
		mapSum = bio1.sum(bio1.axsObservation);
		for(String s:mapSum.keySet()){
			assertEquals(mapSumCorrect.get(s),mapSum.get(s),0.0000000001);
		}
		mapSum = bio1.sum(bio1.axsSample);
		for(String s:mapSum.keySet()){
			assertEquals(mapSumCorrect.get(s),mapSum.get(s),0.0000000001);
		}
	}
	@Test
	public void takeRandomSubset_SubsetIsTaken_TableHasCorrectSize(){
		
		try {
			bio1.takeRandomSubset(3, bio1.axsObservation);
		} catch (Exception e) {
			fail(e.getMessage());
		}
		assertEquals(3,bio1.axsObservation.size());
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
		try {
			bio1.takeRandomSubset(4, bio1.axsObservation);
		} catch (Exception e) {
			fail(e.getMessage());
		}
		assertEquals(4,bio1.axsObservation.size());
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);	
	}
	
	@Test
	public void AxisaddMetadata_MetadataIsAdded_MetadataIsCorrect(){
		
		//mapMetadata = metadata map
		
		HashMap<String,HashMap<String,String>> mapMetadata;
		
		mapMetadata = new HashMap<String,HashMap<String,String>>();
		for(int i=1;i<=5;i++){
			mapMetadata.put("GG_OTU_" + i, new HashMap<String,String>());
			mapMetadata.get("GG_OTU_" + i).put("TestKey",Integer.toString(i));
		}
		bio1.axsObservation.addMetadata(mapMetadata);
		for(String s:mapMetadata.keySet()){
			assertEquals(mapMetadata.get(s).get("TestKey"),bio1.axsObservation.getMetadata(s).get("TestKey"));
		}
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);	
	}
	
	//TODO need to add test for AxisaddMetadataFromTextFile
	
	@Test
	public void AxisgetID_IDGotten_IDCorrect(){
		for(int i=1;i<=5;i++){
			assertEquals("GG_OTU_" + i,bio1.axsObservation.getID(i-1));
		}
	}
	
	@Test
	public void AxisgetIDs_IDsGotten_IDsCorrect(){
		
		//rgsIDsCorrect = correct set of ids
		
		String[] rgsIDsCorrect;
		
		rgsIDsCorrect = new String[]{"GG_OTU_1","GG_OTU_2","GG_OTU_3","GG_OTU_4","GG_OTU_5"};
		assertEquals(rgsIDsCorrect.length,bio1.axsObservation.getIDs().size());
		for(int i=0;i<rgsIDsCorrect.length;i++){
			assertTrue(bio1.axsObservation.getIDs().contains(rgsIDsCorrect[i]));
		}
	}
	
	@Test
	public void AxisgetIndex_IndexGotten_IndexCorrect(){
		for(int i=1;i<=5;i++){
			assertEquals(i-1,bio1.axsObservation.getIndex("GG_OTU_"+i));
		}
	}
	
	//TODO need to write check for AxisgetMetadata
	
	@Test
	public void AxisgetMetadataKeys_KeysGotten_KeysCorrect(){
		
		//rgs1 = correct set of keys
		
		String[] rgs1;
		
		rgs1 = new String[]{"BODY_SITE","BarcodeSequence","Description","LinkerPrimerSequence"};
		assertEquals(rgs1.length,bio1.axsSample.getMetadataKeys().size());
		for(int i=0;i<rgs1.length;i++){
			
			assertTrue(bio1.axsSample.getMetadataKeys().contains(rgs1[i]));
		}
	}
	
	//TODO need to write check for getObjects
	
	@Test
	public void AxishasMetadataField_FieldsExist_ReturnsTrue(){
		
		//rgs1 = set of fields
		
		String[] rgs1;
		
		rgs1 = new String[]{"BODY_SITE","BarcodeSequence","Description","LinkerPrimerSequence"};
		
		for(int i=0;i<rgs1.length;i++){
			assertTrue(bio1.axsSample.hasMetadataField(rgs1[i]));
		}
	}
	
	@Test
	public void AxishasMetadataField_FieldsDoesNotExist_ReturnsFalse(){
		
		//rgs1 = set of fields
		
		String[] rgs1;
		
		rgs1 = new String[]{"BODY_SITE1","BarcodeSequence1","Description1","LinkerPrimerSequence1"};
		
		for(int i=0;i<rgs1.length;i++){
			assertFalse(bio1.axsSample.hasMetadataField(rgs1[i]));
		}
	}
	
	@Test
	public void AxisremoveAllMetadata_MetadataRemoved_NoMetadataRemains(){
		
		bio1.axsSample.removeAllMetadata();
		for(String s:bio1.axsSample.getIDs()){
			assertEquals(0,bio1.axsSample.getMetadata(s).size());
		}
		bio1 = new BiomIO(sTestDataDir + "/" + sTestFile);
	}
	
	@Test
	public void Axissize_SizeChecked_SizeCorrect(){
		assertEquals(5,bio1.axsObservation.size());
	}
	
	private void checkTableIsCorrect(){
		checkAxisIDsAreCorrect(rgsCorrectObservationIDs,bio1.axsObservation);
		checkAxisIDsAreCorrect(rgsCorrectSampleIDs,bio1.axsSample);
		checkDataAreCorrect(rgdCorrectData);
		
		//TODO also check metadata for each observation and sample
	}
	
	private void checkAxisIDsAreCorrect(String[] rgsCorrectIDs, BiomIO.Axis axs1){
		assertEquals(rgsCorrectIDs.length,axs1.size());
		for(int i=0;i<rgsCorrectIDs.length;i++){
			assertTrue(axs1.getIDs().contains(rgsCorrectIDs[i]));
		}
	}
	
	private void checkDataAreCorrect(double[][] rgdCorrectData){
		for(int i=0;i<rgdCorrectData.length;i++){
			for(int j=0;j<rgdCorrectData[0].length;j++){
				assertEquals(rgdCorrectData[i][j],bio1.getValueByIndices(i, j),0.0000001);
			}
		}
	}
	
	private void checkTableMapIsCorrect(HashMap<String,Double> mapCorrectData, HashMap<String,Double> mapBIOMData){
		assertEquals(mapCorrectData.size(),mapBIOMData.size());
		for(String s:mapCorrectData.keySet()){
			assertTrue(mapBIOMData.containsKey(s));
			assertEquals(mapCorrectData.get(s),mapBIOMData.get(s),0.000000001);
		}
	}
}
