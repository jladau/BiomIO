package edu.ucsf.io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

import ucar.ma2.Array;
import ucar.nc2.Attribute;
import ucar.nc2.Group;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

/**
 * Class for biom 2.1 input/output. 
 * @author Joshua Ladau <br/>
 * 		   jladau@gmail.com
 */
public class BiomIO {

	/**Gives the number of non-zero elements in the table.**/
	public int iNNZ;
	
	/**Identifier of the table.**/
	public String sID;
	
	/**Type of table.**/
	public String sType;
	
	/**A string with a static URL providing format details.**/
	public String sFormatURL;
	
	/**The version of the current biom format, major and minor.**/
	public ArrayList<Integer> lstFormatVersion;
	
	/**Package and revision that built the table.**/
	public String sGeneratedBy;
	
	/**Date the table was built (ISO 8601 format).**/
	public String sCreationDate;
	
	/**Sample axis (columns).**/
	public Axis axsSample;

	/**Observation axis (rows). **/
	public Axis axsObservation;
	
	/**Sparse matrix object.**/
	private SparseMatrix spm1;
	

	/**
	 * Constructor that performs specified initial operations on the BIOM table.
	 * @param sBiomPath Absolute path to BIOM file.
	 * @param mapOptions Selected options. Key value pairings can include:
	 *                   <p>
	 *                   <ul>
	 *                   <li>sTaxonRank [string] = Taxonomic units on which to collapse table. Accepted values are "kingdom", "phylum", "class", "order", "family", "genus", "species", or "otu." The value of "otu" will cause table to not be modified.
	 *                   <p>
	 *                   <li>sSampleMetadataPath [string] = Path to text file containing sample metadata formatted according to http://biom-format.org/documentation/adding_metadata.html. For use if BIOM file does not contain metadata. Must include "id" field giving sample IDs.
	 *                   <p>
	 *                   <li>rgsSampleMetadataKeys [list of strings] = Sample metadata keys to load from metadata file. A comma-delimited list.
	 *                   <p>
	 *                   <li>sObservationMetadataPath [string] = Path to text file containing observation metadata formatted according to http://biom-format.org/documentation/adding_metadata.html. For use if BIOM file does not contain metadata. Must include "id" field giving observation IDs.
	 *                   <p>
	 *                   <li>rgsObservationMetadataKeys [array of strings] = Observations metadata keys to load from metadata file. A comma-delimited list.
	 *                   <p>
	 *                   <li>sSamplesToKeepPath [string] = Path to file with list of samples to keep. File should contain a list of sample names.
	 *                   <p>
	 *                   <li>sObservationsToKeepPath [string] = Path to file with list of observations to keep. File should contain a list of observation names.
	 *                   <p>
	 *                   <li>rgsRequiredObservationMetadata [array of strings] = Comma-delimited list of observation metadata keys. Samples lacking data for one or more of these keys will be omitted.
	 *                   <p>
	 *                   <li>rgsRequiredSampleMetadata [list of strings] = Comma-delimited list of sample metadata keys. Samples lacking data for one or more of these keys will be omitted.
	 *                   <p>
	 *                   <li>iRandomSampleSubsetSize [integer] = Number of randomly chosen samples to use. Useful for analyzing large data tables quickly.
	 *                   <p>
	 *                   <li>iRandomObservationSubsetSize [integer] = Number of randomly chosen observations to use. Useful for analyzing large data tables quickly.
	 *                   <p>
	 *                   <li>bCheckRarefied [boolean] = Flag for whether to check for rarefaction. If enabled and table is not rarefied, error will be thrown.
	 *                   <p>
	 *                   <li>bNormalize [boolean] = Flag for whether to normalize within each sample so that entries total to 1.
	 *                   <p>
	 *                   <li>iPrevalenceMinimum [integer] = Minimum prevalence: observations that occur in fewer samples will be omitted from analysis.
	 *                   <p>
	 *                   <li>bPresenceAbsence [boolean] = Flag for whether data should be reduced to presence-absence data.
	 *                   <p>
	 *                   <li>iRarefactionTotal [integer] = Total count to which to rarefy samples.
	 *                   </ul>
	 */
	
	/**
	 * Tests for bug in netcdf library
	 */
	/*
	public void checkStringReadError(){
	
		NetcdfFile fil2 = null;
		
		//loading file
		try{
			fil2 = NetcdfFile.open("/home/jladau/Desktop/AG_100nt_even10k.biom");
		}catch(Exception e){
			e.printStackTrace();
		}
		
		//loading observation/ids variable 
		//Group g1 = fil2.getRootGroup();
		//Group g2 = g1.findGroup("sample");
		//Variable v1 = g2.findVariable("ids");
		Variable v1 = fil2.findVariable("observation/metadata/taxonomy");
		
		//reading string array
		try{
			Array ary1 = v1.read();
			int[] rgiTEMP = ary1.getShape();
			String[][] rgs1 = (String[][]) ary1.copyToNDJavaArray();
			System.out.println(rgs1.length);
			System.out.println(rgs1[0].length);
			System.out.println(rgs1[3][4]);
			System.out.println("HERE");
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	*/
	
	/**
	 * Internal constructor.
	 * @param axsObservation Observation axis.
	 * @param axsSample Sample axis.
	 * @param spm1 Sparse matrix.
	 */
	private BiomIO(Axis axsObservation, Axis axsSample, SparseMatrix spm1){
		this.axsSample = axsSample;
		this.axsObservation = axsObservation;
		this.spm1 = spm1;
		iNNZ = spm1.getNNZ();
	}
	
	/**
	 * Constructor.
	 * @param sBiomPath Absolute path to BIOM file.
	 */
	public BiomIO(String sBiomPath){
		
		//fil1  = File object
		
		NetcdfFile fil1;
		
		//loading file
		fil1 = null;
		try {
			fil1 = NetcdfFile.open(sBiomPath);
			//fil1 = NetcdfDataset.openFile(sBiomPath,null);
		}catch(IOException e){
			e.printStackTrace();
		}
		
		//loading global variables
		loadGlobalVariables(fil1);
		
		//loading axes
		loadAxes(fil1);
		
		//loading taxonomy metadata
		loadTaxonomicMetadata(fil1);
		
		//loading other metadata
		loadMetadata(fil1);
		
		//loading sparse matrix object
		loadSparseMatrix(fil1);
		
		//closing file
		close(fil1);
	}
	
	public BiomIO(String sBiomPath, Map<String,String> mapOptions) throws Exception{
		
		//running constructor
		this(sBiomPath);
		
		//collapsing by taxon if requested
		if(mapOptions.containsKey("sTaxonRank") && axsObservation.hasMetadataField("taxonomy") && !mapOptions.get("sTaxonRank").equals("otu")){
			System.out.println("Collapsing by " + mapOptions.get("sTaxonRank") + "...");
			collapse(mapOptions.get("sTaxonRank"), axsObservation, false);
		
		}
	
		//loading sample metadata from text file if appropriate
		if(mapOptions.containsKey("sSampleMetadataPath")){
			System.out.println("Loading sample metadata from text file...");
			axsSample.removeAllMetadata();
			axsSample.addMetadataFromTextFile(mapOptions.get("sSampleMetadataPath"),mapOptions.get("rgsSampleMetadataKeys").split(","));
		}
		
		//loading observation metadata from text file if appropriate
		if(mapOptions.containsKey("sObservationMetadataPath")){
			System.out.println("Loading observation metadata from text file...");
			axsObservation.removeAllMetadata();
			axsObservation.addMetadataFromTextFile(mapOptions.get("sObservationMetadataPath"),mapOptions.get("rgsObservationMetadataKeys").split(","));
		}
		
		//filtering samples by file
		if(mapOptions.containsKey("sSamplesToKeepPath")){
			System.out.println("Filtering samples by listed file...");
			this.filterFromFile(mapOptions.get("sSamplesToKeepPath"), axsSample);
		}
		
		//rarefying
		if(mapOptions.containsKey("iRarefactionTotal")){
			System.out.println("Rarefying...");
			this.rarefy(Integer.parseInt(mapOptions.get("iRarefactionTotal")));
		}
		
		//checking for rarefaction
		if(mapOptions.containsKey("bCheckRarefied") && Boolean.parseBoolean(mapOptions.get("bCheckRarefied"))){
			System.out.println("Checking whether samples are rarefied...");
			if(checkRarefied()==false){
				System.out.println("Samples are not rarefied. Exiting.");
				throw new Exception();
			}
		}
		
		//filtering observations by file
		if(mapOptions.containsKey("sObservationsToKeepPath")){
			System.out.println("Filtering observations by listed file...");
			this.filterFromFile(mapOptions.get("sObservationsToKeepPath"), axsObservation);
		}

		//filtering observations by no metadata
		if(mapOptions.containsKey("rgsRequiredObservationMetadata")){
			System.out.println("Filtering observations without required metadata...");
			filterByNoMetadata(mapOptions.get("rgsRequiredObservationMetadata").split(","), axsObservation);
		}
		
		//filtering samples by no metadata
		if(mapOptions.containsKey("rgsRequiredSampleMetadata")){
			System.out.println("Filtering samples without required metadata...");
			filterByNoMetadata(mapOptions.get("rgsRequiredSampleMetadata").split(","), axsSample);
		}
		
		//loading random subset of samples
		if(mapOptions.containsKey("iRandomSampleSubsetSize")){
			System.out.println("Taking random subset of samples...");
			takeRandomSubset(Integer.parseInt(mapOptions.get("iRandomSampleSubsetSize")), axsSample);			
		}
		
		//normalizing to relative abundance
		if(mapOptions.containsKey("bNormalize") && Boolean.parseBoolean(mapOptions.get("bNormalize"))){
			System.out.println("Normalizing to relative abundance...");
			normalize();
		}
		
		//filtering by prevalence
		if(mapOptions.containsKey("iPrevalenceMinimum")){
			System.out.println("Filtering observations by prevalence...");
			filterByPrevalence(Integer.parseInt(mapOptions.get("iPrevalenceMinimum")));
		}

		//loading random subset of observations
		if(mapOptions.containsKey("iRandomObservationSubsetSize")){
			System.out.println("Taking random subset of observations...");
			takeRandomSubset(Integer.parseInt(mapOptions.get("iRandomObservationSubsetSize")), axsObservation);
		}
			
		//reducing to presence-absence
		if(mapOptions.containsKey("bPresenceAbsence") && Boolean.parseBoolean(mapOptions.get("bPresenceAbsence"))){
			System.out.println("Converting table to presence-absence data...");
			convertToPresenceAbsence();
		}
	}
	
	/**
	 * Checks whether samples are rarefied.
	 * @return True if all samples have same sum; false otherwise. 
	 */
	public boolean checkRarefied(){
		
		//map1 = map of sums within samples
		//d1 = first sum value
		
		HashMap<String,Double> map1;
		double d1;
		
		map1 = this.sum(axsSample);
		d1 = Double.NaN;
		for(String s:map1.keySet()){
			if(Double.isNaN(d1)){
				d1 = map1.get(s);
			}else{
				if(map1.get(s)!=d1){
					return false;
				}
			}
		}
		return true;
	}
	
	/**
	 * Clears metadata.
	 */
	private void clearMetadata(){
		iNNZ = -9999;
		sID = null;
		sType = null;
		sFormatURL = null;
		lstFormatVersion = null;
		sGeneratedBy = null;
		sCreationDate = null;
	}
	
	/**
	 * Closes netcdf file object.
	 * @param fil1 Netcdf file object.
	 */
	private void close(NetcdfFile fil1){
		if(fil1!=null){ 
			try{
				fil1.close();
			}catch(IOException e){
				e.printStackTrace();
			}
	    }
	}

	
	/**
	 * Resamples samples axis with replacement. Useful for bootstrapping.
	 * @param iRandomSeed Random seed to use.
	 * @return Resampled table. New samples have the same names as original samples with suffixes '.1', '.2', etc.
	 */
	public BiomIO resampleWithReplacement(int iRandomSeed){
		
		//axsNew = new axis
		//spmNew = new sampling matrix
		
		Axis axsNew;
		SparseMatrix spmNew;
		
		axsNew = axsSample.resampleWithReplacement(iRandomSeed);
		spmNew = spm1.resample(axsNew.mapResample);
		return new BiomIO(axsObservation,axsNew,spmNew);
	}
	
	/**
	 * Collapses table (sums entries) by specified metadata field: axis elements are combined if they share a given metadata value.
	 * @param sMetadataKey Key to use for metadata collapsing (elements that share the same value of this metadata will be combined).
	 * @param axs1 Axis to collapse along.
	 * @param bOutputNew Output new BIOM object (versus modify current object).
	 * @return A new BIOM object with axis collapsed by the specified metadata.
	 */
	public BiomIO collapse(String sMetadataKey, Axis axs1, boolean bOutputNew){
		
		//axsSampleNew = new axis
		//spmNew = new sampling matrix
		
		Axis axsNew;
		SparseMatrix spmNew;
		
		if(!axs1.hasMetadataField(sMetadataKey)){
			return this;
		}
		axsNew = axs1.collapse(sMetadataKey);
		if(axsNew.sName.equals("sample")){
			spmNew = spm1.collapse(null, axsNew.mapCollapse);
			if(bOutputNew){	
				return new BiomIO(axsObservation,axsNew,spmNew);
			}else{
				this.axsSample = axsNew;
				this.spm1 = spmNew;
				iNNZ = spmNew.getNNZ();
			}
		}else if(axsNew.sName.equals("observation")){
			spmNew = spm1.collapse(axsNew.mapCollapse, null);
			if(bOutputNew){
				return new BiomIO(axsNew,axsSample,spmNew);
			}else{
				this.axsObservation = axsNew;
				this.spm1 = spmNew;
				iNNZ = spmNew.getNNZ();
			}
		}
		return null;
	}

	/**
	 * Converts table to presence-absence data.
	 */
	public void convertToPresenceAbsence(){
		for(String s:axsObservation.getIDs()){
			for(String t:axsSample.getIDs()){
				if(spm1.getPresenceAbsence(s, t)==1){
					spm1.setValue(s, t, 1);
				}
			}
		}
	}
	
	/**
	 * Checks if BIOM objects are the same.
	 * @param obj1 Object to compare to.
	 * @return True if BIOM objects share the same arrays of data and observation and sample IDs; false otherwise. Note: equality of metadata is not checked.
	 */
	public boolean equals(Object obj1){
		
		//bio1 = BIOM object
		
		BiomIO bio1;
		
		if(obj1 instanceof BiomIO){
			bio1 = (BiomIO) obj1;
		}else{
			return false;
		}
		
		if(axsObservation.size()!=bio1.axsObservation.size()){
			return false;
		}
		if(axsSample.size()!=bio1.axsSample.size()){
			return false;
		}
		for(String s:axsObservation.getIDs()){
			if(!bio1.axsObservation.getIDs().contains(s)){
				return false;
			}
		}
		for(String s:axsSample.getIDs()){
			if(!bio1.axsSample.getIDs().contains(s)){
				return false;
			}
		}
		for(String s:axsObservation.getIDs()){
			for(String t:axsSample.getIDs()){
				if(this.getValueByIDs(s, t)!=bio1.getValueByIDs(s,t)){
					return false;
				}
			}
		}
		return true;
	}
	
	/**
	 * Filter a table based on an iterable.
	 * @param setIDsToKeep Set of IDs of axis elements (i.e., samples or observations) to keep.
	 * @param axs1 Axis from which to remove axis elements (i.e., samples or observations).
	 */
	public void filter(HashSet<String> setIDsToKeep, Axis axs1) throws Exception{
		
		//checking if empty
		if(setIDsToKeep==null || setIDsToKeep.size()==0){
			throw new Exception("Filter will not allow any elements to pass.");
		}
		
		//filtering data matrix	
		if(axs1.sName.equals("sample")){	
			spm1.filterColumns(setIDsToKeep);
		}else if(axs1.sName.equals("observation")){
			spm1.filterRows(setIDsToKeep);
		}
		
		//filtering axes
		axs1.filter(setIDsToKeep);
		
		//updating metadata
		this.clearMetadata();
		
		//updating nozero count
		iNNZ = spm1.getNNZ();
	}

	/**
	 * Filters samples without metadata.
	 * @param rgsMetadataKeys Metadata to check for; if axis element lacks one or more of these keys then the element will be removed.
	 */
	public void filterByNoMetadata(String[] rgsMetadataKeys, Axis axs1) throws Exception{
		
		//setKeep = set of IDs to keep
		
		HashSet<String> setKeep;
		
		setKeep = axs1.getIDs();
		for(AxisObject a:axs1.getObjects()){
			for(int i=0;i<rgsMetadataKeys.length;i++){
				if(!a.hasMetadata(rgsMetadataKeys[i])){
					setKeep.remove(a.sID);
					break;
				}
			}
		}
			
		try{
			this.filter(setKeep, axs1);
		}catch(Exception e){
			System.out.println("No " + axs1.sName + "s with required metadata. Exiting.");
			throw e;
		}
	}
	
	/**
	 * Removes observations that have prevalence below specified level.
	 * @param iPrevalenceMin Minimum prevalence to keep.
	 */
	public void filterByPrevalence(int iPrevalenceMin) throws Exception{
		
		//map1 = returns prevalence of specified observation
		//setKeep = set of observations to keep
		
		HashMap<String,Integer> mapPrevalence;
		HashSet<String> setKeep;
		
		mapPrevalence = this.getNonzeroCounts(this.axsObservation);
		setKeep = new HashSet<String>(mapPrevalence.size());
		for(String s:mapPrevalence.keySet()){
			if(mapPrevalence.get(s)>=iPrevalenceMin){
				setKeep.add(s);
			}
		}
		try{
			this.filter(setKeep, axsObservation);
		}catch(Exception e){
			System.out.println("No observations passed prevalence filter. Exiting.");
			throw e;		
		}
	}
	
	/**
	 * Filters axis based on list of IDs stored in a file.
	 * @param sElementsToKeepPath Path to data file with list of element IDs to keep. File should contain a list of sample or observation IDs.
	 * @param axs1 Axis along which to filter. Either axsSample or axsObservation.
	*/
	public void filterFromFile(String sElementsToKeepPath, Axis axs1) throws Exception{
		
		//setOut = output
		//bfr1 = buffered reader
		//s1 = current line
		
		HashSet<String> setOut;
		BufferedReader bfr1;
		String s1;
		
		setOut = new HashSet<String>(1000);
		try {
			bfr1 = new BufferedReader(new FileReader(sElementsToKeepPath));
			while((s1 = bfr1.readLine()) != null) {
				setOut.add(s1);
			}
			bfr1.close();
		}catch (Exception e){
			e.printStackTrace();
		}
		try{
			filter(setOut,axs1);
		}catch(Exception e){
			System.out.println("No " + axs1.sName + "s in included list. Exiting.");
			throw e;
		}
	}
	
	/**
	 * Gets row or column from matrix.
	 * @param axs1 Axis from which to get row or column.
	 * @param sID ID of observation or sample values to obtain.
	 * @return Map between axis elements (i.e., observation or sample IDs) and values in specified row or column.
	 */
	public HashMap<String,Double> getItem(Axis axs1, String sID){
		
		//mapOut = output
		
		HashMap<String,Double> mapOut;
		
		mapOut = new HashMap<String,Double>(axs1.size());
		if(axs1.sName.equals("sample")){
			for(AxisObject a:axsObservation.getObjects()){
				mapOut.put(a.sID, spm1.getValue(a.sID, sID));
			}
		}else if(axs1.sName.equals("observation")){
			for(AxisObject a:axsSample.getObjects()){
				mapOut.put(a.sID, spm1.getValue(sID, a.sID));
			}
		}
		return mapOut;
	}
	
	/**
	 * Get mean values for an element.
	 * @param axs1 Axis along which to find mean.
	 * @param sID ID of axis element for which to find mean.
	 * @return Mean value for element.
	 */
	public double getMean(Axis axs1, String sID){
		
		//dOut = output
		//dCount = count
		
		double dOut;
		double dCount;
		
		//initializing output
		dOut = 0.;
		dCount = 0.;
		
		//loading values
		if(axs1.sName.equals("sample")){
			for(AxisObject a:axsObservation.getObjects()){
				dOut+=this.getValueByIDs(a.sID, sID);
				dCount++;
			}
		}else if(axs1.sName.equals("observation")){
			for(AxisObject a:axsSample.getObjects()){
				dOut+=this.getValueByIDs(sID, a.sID);
				dCount++;
			}	
		}
		return dOut/dCount;
	}
	
	/**
	 * Get mean values along an axis.
	 * @param axs1 Axis for which to find means.
	 * @return Mean values.
	 */
	public HashMap<String,Double> getMeans(Axis axs1){
		
		//map1 = output
		
		HashMap<String,Double> map1;
		
		//initializing output
		map1 = new HashMap<String,Double>(axs1.size());
		
		//saving results
		for(String s:axs1.getIDs()){
			map1.put(s, this.getMean(axs1,s));
		}
		
		//outputting
		return map1;
	}
	
	/**
	 * Get number of non-zero entries for an axis element.
	 * @param axs1 Axis for which to find non-zero count.
	 * @param sID ID of element for which to find non-zero count.
	 * @return Number of non-zero entries for element.
	 */
	public int getNonzeroCount(Axis axs1, String sID){
		
		//iOut = output
		
		int iOut;
		
		//initializing output
		iOut = 0;
		
		//loading values
		if(axs1.sName.equals("sample")){
			for(AxisObject a:axsObservation.getObjects()){
				if(this.getValueByIDs(a.sID, sID)>0){
					iOut++;
				}
			}
		}else if(axs1.sName.equals("observation")){
			for(AxisObject a:axsSample.getObjects()){
				if(this.getValueByIDs(sID, a.sID)>0){
					iOut++;
				}
			}	
		}
		return iOut;
	}
	
	/**
	 * Get nonzero summaries about an axis.
	 * @param axs1 Axis for which to find counts.
	 * @return Number of non-zero counts.
	 */
	public HashMap<String,Integer> getNonzeroCounts(Axis axs1){
		
		//map1 = output
		
		HashMap<String,Integer> map1;
		
		//initializing output
		map1 = new HashMap<String,Integer>(axs1.size());
		
		//saving results
		for(String s:axs1.getIDs()){
			map1.put(s, this.getNonzeroCount(axs1,s));
		}
		
		//outputting
		return map1;
	}
	
	/**
	 * Gets richness observed in each sample.
	 * @return Richnesses observed in each sample
	 */
	public HashMap<String,Double> getRichness(){
		
		//map1 = output map in integer format
		//map1 = output map in double format
		
		HashMap<String,Integer> map1;
		HashMap<String,Double> map2;
		
		map1 = getNonzeroCounts(axsSample);
		map2 = new HashMap<String,Double>(map1.size());
		for(String s:map1.keySet()){
			map2.put(s, (double) map1.get(s));
		}
		return map2;
	}
	
	/**
	 * Gets Shannon diversity
	 * @return Shannon diversity for each sample
	 */
	public HashMap<String,Double> getShannon(){
		
		//mapOut = output
		//d1 = current shannon diversity
		//d2 = current value
		
		HashMap<String,Double> mapOut;
		double d1;
		double d2;
		
		mapOut = new HashMap<String,Double>(axsSample.size());
		for(String s:axsSample.getIDs()){
			d1 = 0;
			for(String t:axsObservation.getIDs()){
				d2 = this.getValueByIDs(t, s);
				if(d2>0){
					d1 -= d2*Math.log(d2);
				}
			}
			mapOut.put(s, d1);
		}
		return mapOut;
	}
	
	/**
	 * Gets shape of table.
	 * @return Number of rows (observation; entry 0) and columns (samples; entry 1).
	 */
	public int[] getShape(){
		return new int[]{axsObservation.size(),axsSample.size()};
	}

	/**
	 * Gets value for a specified observation and sample.
	 * @param sObservationID ID of observation.
	 * @param sSampleID ID of sample.
	 * @return Value for given observation and sample.
	 */
	public double getValueByIDs(String sObservationID, String sSampleID){
		return spm1.getValue(sObservationID, sSampleID);
	}
	
	/**
	 * Gets value in specified row and column of matrix.
	 * @param iRow Index of observation.
	 * @param iCol Index of sample.
	 * @return Value for given observation and sample.
	 */
	public double getValueByIndices(int iRow, int iCol){
		
		//sRowID = row ID
		//sColID = column ID
				
		String sRowID;
		String sColID;
		
		//loading row and column indices
		sRowID = this.axsObservation.getID(iRow);
		sColID = this.axsSample.getID(iCol);
		
		//returning result
		return getValueByIDs(sRowID,sColID);
	}
	
	
	/**
	 * Loads axes.
	 * @param fil1 Netcdf file object.
	 */
	private void loadAxes(NetcdfFile fil1){
		
		//initializing axes
		axsSample = new Axis("sample", fil1.findVariable("sample/ids"));
		axsObservation = new Axis("observation", fil1.findVariable("observation/ids"));
	}

	/**
	 * Loads top-level attributes.
	 * @param fil1 Netcdf file object.
	 */
	private void loadGlobalVariables(NetcdfFile fil1){
		
		//lst1 = list of global attributes
		
		List<Attribute> lst1;
		lst1 = fil1.getGlobalAttributes();
		for(int i=0;i<lst1.size();i++){
			if(lst1.get(i).getShortName().equals("id")){
				sID=lst1.get(i).getStringValue();
			}
			if(lst1.get(i).getShortName().equals("type")){
				sType=lst1.get(i).getStringValue();
			}
			if(lst1.get(i).getShortName().equals("format-url")){
				sFormatURL=lst1.get(i).getStringValue();
			}
			if(lst1.get(i).getShortName().equals("format-version")){
				lstFormatVersion=new ArrayList<Integer>();
				for(int k=0;k<lst1.get(i).getValues().getSize();k++){
					lstFormatVersion.add(lst1.get(i).getValues().getInt(k));
				}
			}
			if(lst1.get(i).getShortName().equals("generated-by")){
				sGeneratedBy=lst1.get(i).getStringValue();
			}
			if(lst1.get(i).getShortName().equals("creation-date")){
				sCreationDate=lst1.get(i).getStringValue();
			}
			if(lst1.get(i).getShortName().equals("nnz")){
				iNNZ=lst1.get(i).getNumericValue().intValue();
			}
		} 
	}
	
	
	/**
	 * Loads all metadata with the exception of taxonomy metadata.
	 * @param fil1 Netcdf file object.
	 */
	private void loadMetadata(NetcdfFile fil1){
		
		//lstVars = list of metadata variables
		//grp1 = root group
		//grp2 = sample or observation group
		//grp3 = metadata group
		
		List<Variable> lstVars;
		Group grp1;
		Group grp2;
		Group grp3;
		
		grp1 = fil1.getRootGroup();
		grp2 = grp1.findGroup("sample");
		grp3 = grp2.findGroup("metadata");
		lstVars=grp3.getVariables();
		for(int i=0;i<lstVars.size();i++){
			loadMetadataFromVariable(lstVars.get(i),axsSample);
		}
		grp2 = grp1.findGroup("observation");
		grp3 = grp2.findGroup("metadata");
		lstVars=grp3.getVariables();
		for(int i=0;i<lstVars.size();i++){
			if(!lstVars.get(i).getFullName().equals("taxonomy")){
				loadMetadataFromVariable(lstVars.get(i),axsObservation);
			}
		}
	}

	/**
	 * Loads metadata from variable.
	 * @param var1 Variable from which to load metadata.
	 * @param axs1 Axis to add metadata to.
	 */
	private void loadMetadataFromVariable(Variable var1, Axis axs1){
		
		//ary1 = array of metadata
		//rgo1 = current metadata
		//rgo2 = current metadata
		
		Array ary1 = null;
		Object rgo1[][];
		Object rgo2[];
		
		try{
			
			//reading string array
			ary1 = var1.read();
			rgo1 = (Object[][]) ary1.copyToNDJavaArray();

			//saving metadata
			for(int i=0;i<rgo1.length;i++){
				sID=axs1.getID(i);
				for(int j=0;j<rgo1[0].length;j++){
					axs1.setMetadata(sID, var1.getShortName() + "." + j, (String) rgo1[i][j]);
				}
			}
		}catch(Exception e){
			
			try{
				rgo2 = (Object[]) ary1.copyTo1DJavaArray();
				
				//saving metadata
				for(int i=0;i<rgo2.length;i++){
					sID=axs1.getID(i);
					axs1.setMetadata(sID, var1.getShortName(), (String) rgo2[i]);
				}
			}catch(Exception e2){
				e2.printStackTrace();
			}
		}
	}
	
	/**
	 * Loads sparse matrix object.
	 * @param fil1 Netcdf file object.
	 */
	private void loadSparseMatrix(NetcdfFile fil1){
		spm1 = new SparseMatrix(fil1.findVariable("observation/matrix/indices"),fil1.findVariable("observation/matrix/indptr"),fil1.findVariable("observation/matrix/data"), axsObservation, axsSample);
	}
	
	/**
	 * Loads taxonomy metadata.
	 * @param fil1 Netcdf file object.
	 */
	private void loadTaxonomicMetadata(NetcdfFile fil1){
		
		//var1 = variable
		//rgs1 = taxonomy
		//ary1 = array of taxonomy data
		//sID = id
		//rgsClades = clades in order
		//rgsTaxa = current taxonomic assignments
		//rgsAliases = aliases
		//sbl1 = current taxon
		//sPrefix = current prefix
		//mapIndex(sPrefix) = returns index for current prefix
		
		String rgsClades[];
		String rgsTaxa[];
		String rgsAliases[];
		StringBuilder sbl1;
		String rgs1[][] = null;
		Variable var1;
		Array ary1;
		String sID;
		String sPrefix;
		HashMap<String,Integer> mapIndex;
		
		var1 = fil1.findVariable("observation/metadata/taxonomy");
		
		rgsClades = new String[]{"kingdom","phylum","class","order","family","genus","species"};
		rgsAliases = new String[]{"k__","p__","c__","o__","f__","g__","s__"};
		mapIndex = new HashMap<String,Integer>();
		mapIndex.put("k__",0);
		mapIndex.put("p__",1);
		mapIndex.put("c__",2);
		mapIndex.put("o__",3);
		mapIndex.put("f__",4);
		mapIndex.put("g__",5);
		mapIndex.put("s__",6);
		try{
			
			//reading string array
			ary1 = var1.read();
			rgs1 = (String[][]) ary1.copyToNDJavaArray();
			
		}catch(Exception e){
			//System.out.println("Taxonomy metadata not found.");
			return;
		}
			
		//checking if taxonomy data found
		if(rgs1==null){
			System.out.println("Warning: taxonomy metadata not found.");
			return;
		}
		
		//looping through observations
		for(int i=0;i<rgs1.length;i++){
			
			sID=axsObservation.getID(i);
			sbl1 = new StringBuilder();
			rgsTaxa = new String[7];
			
			//loading taxonomy data
			for(int j=0;j<rgs1[0].length;j++){
				
				//cleaning kingdom header
				//if(rgs1[i][j].contains("k__")){
				//	rgs1[i][j]="k__" + rgs1[i][j].split("_")[2];
				//}
				
				//loading taxon
				//if(!rgs1[i][j].equals("") && !rgs1[i][j].toLowerCase().equals("unclassified")){
				//	sPrefix = rgs1[i][j].substring(0, 3);
				//	rgsTaxa[mapIndex.get(sPrefix)]=rgs1[i][j];
				//}else{
				//	rgsTaxa[j]=rgsAliases[j];
				//}
				
				if(!rgs1[i][j].equals("") && !rgs1[i][j].toLowerCase().equals("unclassified")){
					if(rgs1[i][j].length()>3){
						sPrefix = rgs1[i][j].substring(0, 3);
						if(mapIndex.containsKey(sPrefix)){
							rgsTaxa[mapIndex.get(sPrefix)]=rgs1[i][j];
						}
					}
				}
			}
			
			//saving empty taxa
			for(int j=0;j<7;j++){
				if(rgsTaxa[j]==null){
					rgsTaxa[j]=rgsAliases[j];
				}
			}
			
			//saving taxonomy data
			sbl1 = new StringBuilder();
			for(int j=0;j<7;j++){
				if(j>0){
					sbl1.append(";");
				}
				sbl1.append(rgsTaxa[j]);
				if(rgsTaxa[j].length()>3){
					this.axsObservation.setMetadata(sID, rgsClades[j], sbl1.toString());
				}else{
					this.axsObservation.setMetadata(sID, rgsClades[j], "unclassified");
				}
			}
			this.axsObservation.setMetadata(sID, "taxonomy", sbl1.toString());
		}
	}
	
	/**
	 * Normalizes counts within samples: i.e., transforms data to relative abundance.
	 */
	public void normalize(){
		
		//mapSum = sum across samples
		
		HashMap<String,Double> mapSum;
		
		mapSum = this.sum(axsSample);
		for(String s:axsSample.getIDs()){
			if(mapSum.get(s)==0){
				continue;
			}
			for(String t:axsObservation.getIDs()){
				spm1.setValue(t, s, spm1.getValue(t, s)/mapSum.get(s));
			}
		}
	}
	
	/**
	 * Prints metadata for given axis.
	 * @param axs1 Axis for which to print metadata.
	 * @return Comma-delimited table.
	 */
	public ArrayList<String> printMetadata(Axis axs1){
		
		//lstKeys = list of metadata keys
		//lstOut = output
		//sbl1 = current line
		
		ArrayList<String> lstOut;
		ArrayList<String> lstKeys;
		StringBuilder sbl1;
		
		//initializing output
		lstOut = new ArrayList<String>(axs1.size()+1);
		
		//loading list of keys
		lstKeys = new ArrayList<String>();
		for(String s:axs1.getMetadataKeys()){
			lstKeys.add(s);
		}
		
		//outputting headers
		sbl1 = new StringBuilder();
		sbl1.append(axs1.sName);
		for(int i=0;i<lstKeys.size();i++){		
			sbl1.append("," + lstKeys.get(i));
		}
		lstOut.add(sbl1.toString());
		
		//outputting data
		for(AxisObject a:axs1.getObjects()){
			sbl1 = new StringBuilder();
			sbl1.append(a.sID);
			for(int i=0;i<lstKeys.size();i++){
				sbl1.append("," + a.getMetadata(lstKeys.get(i)));
			}
			lstOut.add(sbl1.toString());
		}
		
		//returning output
		return lstOut;
	}
	
	/**
	 * Prints table.
	 * @return Comma-delimited table.
	 */
	public ArrayList<String> printTable(){
		
		//lstOut = output
		//sbl1 = current line
		
		ArrayList<String> lstOut;
		StringBuilder sbl1;
		
		//outputting headers
		lstOut = new ArrayList<String>(axsObservation.size()+2);
		lstOut.add("# Constructed from biom file");
		sbl1 = new StringBuilder();
		sbl1.append("#OTU ID");
		for(int j=0;j<axsSample.size();j++){
			sbl1.append("," + axsSample.getID(j));
		}
		lstOut.add(sbl1.toString());
		
		//outputting data
		for(int i=0;i<axsObservation.size();i++){
			sbl1 = new StringBuilder();
			sbl1.append(axsObservation.getID(i));
			for(int j=0;j<axsSample.size();j++){
				sbl1.append("," + this.getValueByIndices(i, j));
			}
			lstOut.add(sbl1.toString());
		}
		
		//returning output
		return lstOut;
	}
	
	
	/**
	 * Rarefies samples to specified total.
	 * @param iTotal Total to rarefy to.
	 */
	public void rarefy(int iTotal) throws Exception{
		this.subsample(iTotal, axsSample);
	}
	
	/**
	 * Rarefies vector.
	 * @param iTotal Total to which to rarefy.
	 * @param iVectorTotal Total observations in vector.
	 * @param mapVector Map in which keys are cumulative counts, values are axis element (sample or observation) IDs.
	 * @return Map in which keys are axis element IDs (i.e., sample or observation IDs), values are frequencies.
	 */
	private HashMap<String,Integer> rarefyVector(int iTotal, TreeMap<Integer,String> mapVector, int iVectorTotal){
		
		//i2 = current count
		//set1 = contains indices for selected values
		//mapOut = output; gives the frequency of each object
		//s1 = current object
		
		int i2;
		HashSet<Integer> set1;
		HashMap<String,Integer> mapOut;
		String s1;
		
		//checking if total for vector is too small
		if(iVectorTotal<iTotal){
			return null;
		}
		
		//getting rarefied indices
		set1 = sampleWithoutReplacement(iTotal,iVectorTotal);
		
		//initializing output
		mapOut = new HashMap<String,Integer>(mapVector.size());
		for(Integer i:set1){
			s1 = mapVector.get(mapVector.floorKey(i));
			if(!mapOut.containsKey(s1)){
				mapOut.put(s1, 1);
			}else{
				i2 = mapOut.get(s1);
				i2++;
				mapOut.put(s1, i2);
			}
		}
		
		//returning result
		return mapOut;
	}
	
	/**
	 * Samples without replacement from the set 0,1,...,n.
	 * @param iSampleSize Total number of elements to sample.
	 * @param iPopulationSize Total number of elements from which to sample.
	 * @return A set of integers giving the subsample.
	 */
	private HashSet<Integer> sampleWithoutReplacement(int iSampleSize, int iPopulationSize){
		
		//setOut = output
		//n, N = Knuth's variable names
		//t = total input records dealt with
		//m = number of items selected so far
		//u = random number 
		
		HashSet<Integer> setOut;
		int n;
	    int N;
	    int t = 0; 
	    int m = 0; 
	    double u;

	    setOut = new HashSet<Integer>(iSampleSize);
	    n = iSampleSize;
	    N = iPopulationSize;
	    
	    while (m < n){
	        u = Math.random(); // call a uniform(0,1) random number generator
	        if( (N - t)*u >= n - m ){
	            t++;
	        }else{
	            setOut.add(t);
	            t++; 
	            m++;
	        }
	    }
		return setOut;
	}

	/**
	 * Randomly subsample (rarefy) without replacement.
	 * @param iTotal Total for each subsampling.
	 * @param axs1 Axis along which to subsample.
	*/
	public void subsample(int iTotal, Axis axs1) throws Exception{
		
		//map1 = current vector of data
		//map2 = current cumulative vector of data
		//map3 = rarefied vector
		//i1 = current cumulant
		//setKeep = set of rows or columns to keep
		//dValue = value being input

		TreeMap<Integer,String> map2;
		HashMap<String,Double> map1;
		HashMap<String,Integer> map3;
		int i1;
		HashSet<String> setKeep;
		double dValue;
		
		setKeep = new HashSet<String>(axs1.size());
		for(AxisObject a:axs1.getObjects()){
			
			//loading vectors
			map1 = this.getItem(axs1, a.sID);
			map2 = new TreeMap<Integer,String>();
			i1=0;
			for(String s:map1.keySet()){
				if(map1.get(s)>0){
					map2.put(i1, s);
					i1+=map1.get(s);
				}
			}
			
			//rarefying vector
			map3 = rarefyVector(iTotal,map2,i1);
			
			//rarefaction not possible
			if(map3!=null){
				
				//updating sparse matrix
				for(String s:map1.keySet()){
					
					if(map3.containsKey(s)){
						dValue=map3.get(s);
					}else{
						dValue=0;
					}
					if(axs1.sName.equals("sample")){
						this.spm1.setValue(s, a.sID, dValue);
					}else if(axs1.sName.equals("observation")){
						this.spm1.setValue(a.sID, s, dValue);
					}
				}
				
				//saving element to keep
				setKeep.add(a.sID);
			}
		}
		
		//filtering matrix to remove rows or column that didn't have enough observations
		try{	
			this.filter(setKeep, axs1);
		}catch(Exception e){
			System.out.println("No " + axs1.sName + "s with enough data for subsampling depth. Exiting.");
			throw e;
		}
	}
	
	/**
	 * Returns row or column sums.
	 * @param axs1 Axis for which to find sums.
	 * @return Sums of observations.
	 */
	public HashMap<String,Double> sum(Axis axs1){
		
		//mapOut = output
		
		HashMap<String,Double> mapOut;
	
		mapOut = new HashMap<String,Double>();
		for(String s:axs1.getIDs()){
			mapOut.put(s, spm1.getMarginalSum(axs1.sName, s));
		}
		return mapOut;
	}

	/**
	 * Keeps random subset of elements from axis. Other elements are removed.
	 * @param iSubsetSize Number of elements to keep.
	 * @param axs1 Axis to consider.
	 */
	public void takeRandomSubset(int iSubsetSize, Axis axs1) throws Exception{
		
		//lst1 = axis ids in randomized order
		//setKeep = set of elements to keep
		
		ArrayList<String> lst1;
		HashSet<String> setKeep;
		
		if(iSubsetSize>axs1.size()){
			return;
		}
		
		if(iSubsetSize>axs1.size()){
			System.out.println("Cannot take subset of size " + iSubsetSize + " from " + axs1.sName + " axis, which has " + axs1.size() + " elments. Exiting.");
			throw new Exception();
		}
		
		lst1 = new ArrayList<String>(axs1.size());
		for(int i=0;i<axs1.size();i++){
			if(spm1.getMarginalSum(axs1.sName, axs1.getID(i))>0){
				lst1.add(axs1.getID(i));
			}
		}
		Collections.shuffle(lst1);
		setKeep = new HashSet<String>(iSubsetSize);
		for(int i=0;i<iSubsetSize;i++){
			setKeep.add(lst1.get(i));
		}
		this.filter(setKeep, axs1);
	}

	/**
	 * Axis object.
	 */
	public class Axis{
		
		//sName = name of axis; either 'sample' or 'observation'
		//lstObjects(i) = returns the ith object
		//mapIndex(sID) = returns the index of the specified object id
		//mapCollapse(iIndexOld) = returns new index (for axes that have been collapsed; null otherwise)
		//setMetadataKeys = metadata keys
		
		/**Name of axis; either "sample" or "observation".**/
		public String sName;
		
		/**Set of all metadata keys.**/
		private HashSet<String> setMetadataKeys;
		
		/**Returns the index of the specified object ID.**/
		private HashMap<String,Integer> mapIndex;
		
		/**Returns new ID (for axes that have been collapsed; null otherwise).**/
		private HashMap<String,String> mapCollapse = null;
		
		/**Returns old ID (for axes that have been resampled; null otherwise).**/
		private HashMap<String,String> mapResample = null;
		
		/**Returns the ith object.**/
		private ArrayList<AxisObject> lstObjects;
	
		/**
		 * Internal constructor.
		 * @param sName Axis name: either "observation" or "sample".
		 * @param mapIndex Index map.
		 * @param lstObjects Objects list.
		 * @param setMetadataKeys Metadata keys.
		 * @param mapCollapse Collapse map.
		 */
		private Axis(String sName, HashMap<String,Integer> mapIndex, ArrayList<AxisObject> lstObjects, HashSet<String> setMetadataKeys, HashMap<String,String> mapCollapse, HashMap<String,String> mapResample){
			this.sName=sName;
			this.mapIndex = mapIndex;
			this.lstObjects = lstObjects;
			this.setMetadataKeys = setMetadataKeys;
			this.mapCollapse = mapCollapse;
			this.mapResample = mapResample;
		}
	
		/**
		 * Constructor.
		 * @param sName Axis name: either "observation" or "sample".
		 * @param var1 Netcdf variable with axis element IDs.
		 */
		private Axis(String sName, Variable var1){
			this.sName=sName;
			this.initializeObjects(var1);
			this.setMetadataKeys = new HashSet<String>();
		}
		
		/**
		 * Take a dictionary of metadata and add it to axis.
		 * @param mapMetadata Map in which keys are axis element IDs, values are maps between metadata headings and values.
		 */
		public void addMetadata(HashMap<String,HashMap<String,String>> mapMetadata){
			for(String s:mapMetadata.keySet()){
				for(String t:mapMetadata.get(s).keySet()){
					this.setMetadata(s, t, mapMetadata.get(s).get(t));
				}
			}
		}
		
		/**
		 * Adds metadata from a text file.
		 * @param sMetadataPath Path to file with metadata. Formatted according to http://biom-format.org/documentation/adding_metadata.html. Must include "id" field.
		 * @param rgsFields Metadata keys to add. If Keys are not found in the metadata file they will not be added.
		 */
		public void addMetadataFromTextFile(String sMetadataPath, String[] rgsFields){
			
			//map1 = map of values
			//bfr1 = buffered reader
			//s1 = current line
			//rgs1 = current line in split format
			//mapCol(sMetadataHeader) = returns the column for given metadata header
			
			HashMap<String,HashMap<String,String>> map1;
			HashMap<String,Integer> mapCol=null;
			BufferedReader bfr1;
			String s1;
			String rgs1[];
			
			map1 = new HashMap<String,HashMap<String,String>>();
			try{
				
				//reading input
				bfr1 = new BufferedReader(new FileReader(sMetadataPath));
				while((s1 = bfr1.readLine()) != null){
					if(!s1.startsWith("#")){
						rgs1=s1.split("\t");					
						if(mapCol==null){
							mapCol = new HashMap<String,Integer>();
							for(int i=0;i<rgs1.length;i++){
								mapCol.put(rgs1[i].toLowerCase(), i);
							}
						}else{
							map1.put(rgs1[mapCol.get("id")], new HashMap<String,String>());
							for(int i=0;i<rgsFields.length;i++){
								if(mapCol.containsKey(rgsFields[i])){
									map1.get(rgs1[mapCol.get("id")]).put(rgsFields[i], rgs1[mapCol.get(rgsFields[i])]);
								}
							}
						}
					}
			    }
				bfr1.close();
			}catch (Exception e){
				e.printStackTrace();
			}
			this.addMetadata(map1);
		}
		
		/**
		 * Collapses axis on specified metadata field.
		 * @param sMetadataKey Metadata field to use for collapse.
		 * @return New axis object with elements combined by metadata key.
		 */
		private Axis collapse(String sMetadataKey){
			
			//lstObjectsNew = new set of objects
			//mapIndexNew = new index map
			//mapCollapseNew = new collapse map
			//sIDNew = new id
			//iCounter = new counter
			//setMetadataKeysNew = new set of metadata keys
			
			int iCounter;
			HashMap<String,Integer> mapIndexNew;
			HashMap<String,String> mapCollapseNew;
			ArrayList<AxisObject> lstObjectsNew;
			String sIDNew;
			HashSet<String> setMetadataKeysNew;
			
			lstObjectsNew = new ArrayList<AxisObject>(lstObjects.size());
			mapIndexNew = new HashMap<String,Integer>();
			mapCollapseNew = new HashMap<String,String>();
			setMetadataKeysNew = new HashSet<String>();
			setMetadataKeysNew.add(sMetadataKey);
			iCounter = 0;
			for(AxisObject a:lstObjects){
				sIDNew = a.getMetadata(sMetadataKey);
				if(!mapIndexNew.containsKey(sIDNew)){
					lstObjectsNew.add(new AxisObject(sIDNew));
					mapIndexNew.put(sIDNew, iCounter);
					iCounter++;
				}
				mapCollapseNew.put(a.sID, sIDNew);
			}
			return new Axis(sName, mapIndexNew, lstObjectsNew, setMetadataKeysNew, mapCollapseNew, null);
		}
		
		/**
		 * Resamples axis with replacement. Useful for bootstrap resampling
		 * @param iRandomSeed Random seed.
		 * @return New axis object with resampled elements.
		 */
		private Axis resampleWithReplacement(int iRandomSeed){
			
			//lstObjectsNew = new set of objects
			//mapIndexNew = new index map
			//mapResampleNew = resampled map
			//sIDNew = new id
			//rnd1 = random number generator
			//i1 = current random index
			//i2 = current counter
			//mapIndex(sID) = returns the current index for string
			//sID = current object ID
			//axo = current new axis object
			
			HashMap<String,Integer> mapIndex;
			int i2;
			int i1;
			Random rnd1;
			HashMap<String,Integer> mapIndexNew;
			HashMap<String,String> mapResampleNew;
			ArrayList<AxisObject> lstObjectsNew;
			String sID;
			AxisObject axo1;
			
			lstObjectsNew = new ArrayList<AxisObject>(lstObjects.size());
			mapIndexNew = new HashMap<String,Integer>();
			rnd1 = new Random(iRandomSeed);
			mapIndex = new HashMap<String,Integer>();
			mapResampleNew = new HashMap<String,String>();
			for(int i=0;i<this.size();i++){
				
				i1 = rnd1.nextInt(this.size());
				sID = lstObjects.get(i1).sID;
				if(!mapIndex.containsKey(sID)){
					mapIndex.put(sID, 1);
				}
				i2 = mapIndex.get(sID);
				i2++;
				mapIndex.put(sID, i2);
				sID = sID + "." + (mapIndex.get(sID)-1);
				axo1 = new AxisObject(sID);
				for(String s:lstObjects.get(i1).getAllMetadata().keySet()){
					axo1.addMetadata(s, lstObjects.get(i1).getAllMetadata().get(s));
				}
				lstObjectsNew.add(axo1);
				mapResampleNew.put(sID, lstObjects.get(i1).sID);
				mapIndexNew.put(sID, i);
			}
			return new Axis(sName, mapIndexNew, lstObjectsNew, setMetadataKeys, null, mapResampleNew);
		}

		/**
		 * Filters axis elements.
		 * @param setIDsToKeep Set of IDs of axis elements to keep.
		 * @throws Exception If no axis elements pass filter.
		 */
		private void filter(HashSet<String> setIDsToKeep) throws Exception{
			
			HashMap<String,Integer> mapIndex;
			ArrayList<AxisObject> lstObjects;
			
			mapIndex = new HashMap<String,Integer>(this.mapIndex.size());
			lstObjects = new ArrayList<AxisObject>(this.lstObjects.size());
			
			for(int i=0;i<this.lstObjects.size();i++){			
				if(setIDsToKeep.contains(this.lstObjects.get(i).sID)){
					lstObjects.add(this.lstObjects.get(i));
					mapIndex.put(this.lstObjects.get(i).sID, lstObjects.size()-1);
				}
			}
			this.mapCollapse=null;
			this.lstObjects = lstObjects;
			this.mapIndex=mapIndex;
			if(lstObjects.size()==0){
				throw new Exception();
			}
		}
		
		
		/**
		 * Gets the ID of the axis element associated with the specified index.
		 * @param iIndex Index of axis element.
		 * @return ID of observation or sample.
		 */
		public String getID(int iIndex){
			return lstObjects.get(iIndex).sID;
		}
		
		/**
		 * Return the ids along the given axis.
		 * @return Set of IDs of elements on axis.
		 */
		public HashSet<String> getIDs(){
			
			//lstOut = output
			
			HashSet<String> setOut;
			
			setOut = new HashSet<String>(this.size());
			for(int i=0;i<this.size();i++){
				setOut.add(this.getID(i));
			}
			return setOut;
		}
		
		/**
		 * Get the index of a specified sample or observation along axis.
		 * @param sID Identity of the sample or observation whose index will be returned.
		 * @return Index of sample or observation.
		 */
		public int getIndex(String sID){
			return mapIndex.get(sID);
		}
		
		
		/**
		 * Get the metadata of the identified sample or observation.
		 * @param iIndex Index of axis element (i.e., sample or observation).
		 * @return Map in which keys are metadata headings, values are metadata values.
		 */
		public HashMap<String,String> getMetadata(int iIndex){
			return lstObjects.get(iIndex).getAllMetadata();
		}
		
		/**
		 * Get the metadata of the identified axis sample or observation.
		 * @param sID ID of axis element (i.e., sample or observation).
		 * @return Map in which keys are metadata headings, values are metadata values.
		 */
		public HashMap<String,String> getMetadata(String sID){
			return lstObjects.get(mapIndex.get(sID)).getAllMetadata();
		}
		
		/**
		 * Returns all metadata keys for axis.
		 * @return Set of all metadata keys for axis.
		 */
		public HashSet<String> getMetadataKeys(){
			return this.setMetadataKeys;
		}
		
		/**
		 * Gets list of axis objects.
		 * @return List of all axis elements.
		 */
		private ArrayList<AxisObject> getObjects(){
			return this.lstObjects;
		}
		
		
		/**
		 * Checks whether elements of axis have specified type of metadata.
		 * @param sKey Type of metadata.
		 * @return True if elements have metadata field; false otherwise.
		 */
		public boolean hasMetadataField(String sKey){
			if(setMetadataKeys.contains(sKey)){
				return true;
			}else{
				return false;
			}
		}
		
		/**
		 * Initializes axis elements.
		 * @param var1 Netcdf variable with IDs of axis elements.
		 */
		private void initializeObjects(Variable var1){
			
			//ary1 = array of data
			//rgs1 = current array of axis object names
			
			Array ary1;
			String rgs1[];
			
			//initializing maps
			mapIndex = new HashMap<String,Integer>(1000);
			lstObjects = new ArrayList<AxisObject>(1000);
	
			//loading axis objects
			try{
				ary1 = var1.read();
				rgs1 = (String[]) ary1.copyToNDJavaArray();
				for(int i=0;i<rgs1.length;i++){
					mapIndex.put(rgs1[i], i);
					lstObjects.add(new AxisObject(rgs1[i]));
				}
			}catch(Exception e){
				e.printStackTrace();
			}
		}
		
		/**
		 * Clears all metadata.
		 */
		public void removeAllMetadata(){
			for(AxisObject a:this.getObjects()){
				a.removeMetadata();
			}
		}
		
		/**
		 * Sets metadata.
		 * @param sID ID of axis element.
		 * @param sKey Metadata field.
		 * @param sValue Metadata value.
		 */
		private void setMetadata(String sID, String sKey, String sValue){
			
			if(!mapIndex.containsKey(sID)){
				return;
			}else{		
				if(lstObjects.get(mapIndex.get(sID)).addMetadata(sKey, sValue)==1){
					setMetadataKeys.add(sKey);
				}
			}
		}
		
		/**
		 * Return the length of an axis.
		 * @return Number of elements along axis.
		 */		
		public int size(){
			return lstObjects.size();
		}
	}

	/**
	 * Element of axis.
	 */
	private class AxisObject{
		
		/**Missing data value.**/
		private static final String MISSING_DATA_VALUE="NA";
		
		/**Metadata map.**/
		private HashMap<String,String> mapMetadata;
		
		/**ID of axis object.**/
		private String sID;
		
		/**
		 * Constructor.
		 * @param sID ID of axis element.
		 */
		private AxisObject(String sID){
			this.sID = sID;
			mapMetadata = new HashMap<String,String>();
		}
		
		/**
		 * Adds metadata value.
		 * @param sKey Metadata field.
		 * @param sValue Metadata value.
		 * @return 1 if successful, 0 otherwise.
		 */
		private int addMetadata(String sKey, String sValue){
			if(sValue!=null && !sValue.equals("") && !sValue.equals(MISSING_DATA_VALUE)){
				mapMetadata.put(sKey, sValue);
				return 1;
			}
			return 0;
		}
		
		/**
		 * Gets all metadata.
		 * @return Map between all metadata fields and values.
		 */
		private HashMap<String,String> getAllMetadata(){
			return mapMetadata;
		}
		
		/**
		 * Gets metadata value.
		 * @param sKey Metadata field.
		 * @return Metadata value for specified field; null if field not found.
		 */
		private String getMetadata(String sKey){
			if(!mapMetadata.containsKey(sKey)){
				return null;
			}else{
				return mapMetadata.get(sKey);
			}
		}
		
		/**
		 * Checks if axis element has metadata with given key.
		 * @param sKey Key.
		 * @return True if metadata key found; false otherwise.
		 */
		private boolean hasMetadata(String sKey){
			if(!mapMetadata.containsKey(sKey)){
				return false;
			}else{
				return true;
			}
		}
		
		/**
		 * Clears metadata.
		 */
		private void removeMetadata(){
			mapMetadata = new HashMap<String,String>();
		}
		
		public String toString(){
			return sID;
		}
	}

	/**
	 * Sparse matrix: used for looking up values.
	 */
	private class SparseMatrix{
		
		/**Returns hashmap: keys are column IDs with occurrences, values data values**/
		private HashMap<String,HashMap<String,Double>> mapValue;
		
		/**Returns row sum for given observation**/
		private HashMap<String,Double> mapRowSum;
		
		/**Returns column sum for given sample**/
		private HashMap<String,Double> mapColSum;
	
		/**
		 * Internal constructor.
		 * @param mapValue Value map for new sparse matrix object.
		 */
		private SparseMatrix(HashMap<String,HashMap<String,Double>> mapValue){
			this.mapValue = mapValue;
			mapColSum=null;
			mapRowSum=null;
		}
		
		/**
		 * Constructor.
		 * @param varColIndices Netcdf variable with column indices.
		 * @param varRowPtrs Netcdf variable with row pointers.
		 * @param varData Netcdf variable with data.
		 * @param axsObservation Observation axis.
		 * @param axsSample Sample axis.
		 */		
		private SparseMatrix(Variable varColIndices, Variable varRowPtrs, Variable varData, Axis axsObservation, Axis axsSample){
			
			//aryColIndices = column indices array
			//aryRowPtrs = row pointers array
			//mapPtr(iIndex) = returns specified pointer
			//sRowID = current row
			//sColID = current column
			//aryData = data
			
			TreeMap<Integer,Integer> mapPtr;
			String sRowID; String sColID;
			Array aryColIndices = null;
			Array aryRowPtrs = null;
			Array aryData = null;
			
			//loading arrays
			try{
				aryColIndices = varColIndices.read();
				aryRowPtrs = varRowPtrs.read();
				aryData = varData.read();
			}catch(Exception e){
				e.printStackTrace();
			}
			
			//loading pointer map
			mapPtr = new TreeMap<Integer,Integer>();
			for(int i=0;i<aryRowPtrs.getSize();i++){
				mapPtr.put(aryRowPtrs.getInt(i),i);
			}
		
			//loading value map
			mapValue = new HashMap<String,HashMap<String,Double>>();
			for(int i=0;i<aryColIndices.getSize();i++){
				sRowID = axsObservation.getID(mapPtr.get(mapPtr.floorKey(i)));
				sColID = axsSample.getID(aryColIndices.getInt(i));
				if(!mapValue.containsKey(sRowID)){
					mapValue.put(sRowID, new HashMap<String,Double>());
				}
				mapValue.get(sRowID).put(sColID, aryData.getDouble(i));
			}
			
			//clearning marginal sums
			mapRowSum=null;
			mapColSum=null;
		}
		
		/**
		 * Collapses matrix.
		 * @param mapRow Map from old observation IDs to new observation IDs; null if no collapsing to be done on observations.
		 * @param mapCol Map from old sample IDs to new sample IDs; null if no collapsing to be done on samples.
		 * @return Sparse matrix object with collapsed rows and columns.
		 */
		private SparseMatrix collapse(HashMap<String,String> mapRow, HashMap<String,String> mapCol){
			
			//mapOut = output value map
			//sRowIDNew = new row
			//sColIDNew = new column
			//d1 = current value
			
			HashMap<String,HashMap<String,Double>> mapOut;
			String sRowIDNew;
			String sColIDNew;
			double d1;
			
			mapOut = new HashMap<String,HashMap<String,Double>>();
			for(String s:mapValue.keySet()){
				if(mapRow!=null){	
					sRowIDNew = mapRow.get(s);
				}else{
					sRowIDNew=s;
				}
				if(!mapOut.containsKey(sRowIDNew)){
					mapOut.put(sRowIDNew, new HashMap<String,Double>());
				}
				for(String t:mapValue.get(s).keySet()){
					if(mapCol!=null){
						sColIDNew = mapCol.get(t);
					}else{
						sColIDNew=t;
					}
					if(mapOut.get(sRowIDNew).containsKey(sColIDNew)){
						d1 = mapOut.get(sRowIDNew).get(sColIDNew);
					}else{
						d1=0;
					}
					d1+=mapValue.get(s).get(t);
					mapOut.get(sRowIDNew).put(sColIDNew,d1);
				}
			}
			return new SparseMatrix(mapOut);
		}
		
		/**
		 * Resamples matrix.
		 * @param mapResample Map from new sample IDs to old sample IDs.
		 * @return Sparse matrix object with resampled columns.
		 */
		private SparseMatrix resample(HashMap<String,String> mapResample){
			
			//mapOut = output value map
			
			HashMap<String,HashMap<String,Double>> mapOut;
			
			mapOut = new HashMap<String,HashMap<String,Double>>();
			for(String s:mapValue.keySet()){
				mapOut.put(s, new HashMap<String,Double>());
				for(String t:mapResample.keySet()){
					if(mapValue.get(s).containsKey(mapResample.get(t))){
						mapOut.get(s).put(t, mapValue.get(s).get(mapResample.get(t)));
					}
				}
			}
			return new SparseMatrix(mapOut);
		}
		
		/**
		 * Filters columns of matrix.
		 * @param setColumnsToKeep Set of IDs of samples to keep.
		 */
		private void filterColumns(HashSet<String> setColumnsToKeep){
			
			//map1 = replacement map
			
			HashMap<String,HashMap<String,Double>> map1;
	
			map1 = new HashMap<String,HashMap<String,Double>>(mapValue.size());
			for(String s:mapValue.keySet()){
				for(String t:mapValue.get(s).keySet()){
					if(setColumnsToKeep.contains(t)){
						if(!map1.containsKey(s)){
							map1.put(s, new HashMap<String,Double>(mapValue.get(s).size()));
						}
						map1.get(s).put(t, mapValue.get(s).get(t));
					}
				}
			}
			mapValue=map1;
			mapColSum=null;
			mapRowSum=null;
		}
		
		/**
		 * Filters rows of matrix.
		 * @param setRowsToKeep Set of IDs of observations to keep.
		 */
		private void filterRows(HashSet<String> setRowsToKeep){
			
			//map1 = replacement map
			
			HashMap<String,HashMap<String,Double>> map1;
	
			map1 = new HashMap<String,HashMap<String,Double>>(setRowsToKeep.size());
			for(String s:setRowsToKeep){
				if(mapValue.containsKey(s)){
					map1.put(s, mapValue.get(s));
				}
			}
			mapValue=map1;
			mapRowSum=null;
			mapColSum=null;
		}
		
		/**
		 * Gets marginal sum.
		 * @param sAxisName Axis for which to get marginal sum.
		 * @param sID Axis element for which to get marginal sum.
		 * @return Marginal sum for axis element.
		 */
		private double getMarginalSum(String sAxisName, String sID){
			if(mapRowSum==null || mapColSum==null){
				loadMarginalSums();
			}
			if(sAxisName.equals("sample")){
				if(mapColSum.containsKey(sID)){
					return mapColSum.get(sID);
				}else{
					return 0.;
				}
			}else if(sAxisName.equals("observation")){
				if(mapRowSum.containsKey(sID)){
					return mapRowSum.get(sID);
				}else{
					return 0.;
				}
			}else{
				return Double.NaN;
			}
		}
		
		/**
		 * Gets number of non-zero values.
		 * @return Number of non-zero values.
		 */
		private int getNNZ(){
			int iOut;
			iOut=0;
			for(String s:mapValue.keySet()){
				iOut+=mapValue.get(s).size();
			}
			return iOut;
		}
		
		/**
		 * Gets presence or absence value.
		 * @param sRowID Observation ID.
		 * @param sColID Sample ID.
		 * @return 1 for presence, 0 for absence.
		 */
		private int getPresenceAbsence(String sRowID, String sColID){
			if(mapValue.containsKey(sRowID)){
				if(mapValue.get(sRowID).containsKey(sColID)){
					return 1;
				}else{
					return 0;
				}
			}else{
				return 0;
			}
		}
	
		/**
		 * Returns value at specified location.
		 * @param sRowID Observation ID.
		 * @param sColID Sample ID.
		 * @return Value in matrix.
		 */
		private double getValue(String sRowID, String sColID){
			if(getPresenceAbsence(sRowID,sColID)==0){
				return 0;
			}else{
				return mapValue.get(sRowID).get(sColID);
			}
		}
		
		/**
		 * Loads marginal sum vectors.
		 */
		private void loadMarginalSums(){
			
			//d1 = current value
			
			double d1;
			
			mapRowSum = new HashMap<String,Double>(mapValue.size());
			mapColSum = new HashMap<String,Double>(mapValue.size());
			
			for(String sRowID:mapValue.keySet()){
				for(String sColID:mapValue.get(sRowID).keySet()){
					
					//updating row sum
					if(!mapRowSum.containsKey(sRowID)){
						mapRowSum.put(sRowID,0.);
					}
					d1=mapRowSum.get(sRowID);
					d1+=mapValue.get(sRowID).get(sColID);
					mapRowSum.put(sRowID, d1);
					
					//updating column sum
					if(!mapColSum.containsKey(sColID)){
						mapColSum.put(sColID, 0.);
					}
					d1=mapColSum.get(sColID);
					d1+=mapValue.get(sRowID).get(sColID);
					mapColSum.put(sColID, d1);
				}
			}
		}
		
		/**
		 * Sets value at specified location.
		 * @param sRowID Observation ID.
		 * @param sColID Sample ID.
		 * @param Value in matrix.
		 */
		private void setValue(String sRowID, String sColID, double dValue){
			
			//initializing sums if appropriate
			if(mapRowSum==null || mapColSum==null){
				loadMarginalSums();
			}
			
			//updating row sum
			if(!mapRowSum.containsKey(sRowID)){
				mapRowSum.put(sRowID, dValue);
			}else{
				mapRowSum.put(sRowID, mapRowSum.get(sRowID)-getValue(sRowID,sColID) + dValue);
			}
			if(mapRowSum.get(sRowID)==0){
				mapRowSum.remove(sRowID);
			}
			
			//updating column sum
			if(!mapColSum.containsKey(sColID)){
				mapColSum.put(sColID, dValue);
			}else{
				mapColSum.put(sColID, mapColSum.get(sColID)-getValue(sRowID,sColID) + dValue);
			}
			if(mapColSum.get(sColID)==0){
				mapColSum.remove(sColID);
			}
			
			if(dValue==0){
				if(this.getPresenceAbsence(sRowID, sColID)==1){
					mapValue.get(sRowID).remove(sColID);
					if(mapValue.get(sRowID).size()==0){
						mapValue.remove(sRowID);
					}	
				}
			}else{
				if(!mapValue.containsKey(sRowID)){
					mapValue.put(sRowID, new HashMap<String,Double>());
				}
				mapValue.get(sRowID).put(sColID, dValue);
			}
		}
	}
}
