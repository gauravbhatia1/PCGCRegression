package estimation;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.stat.StatUtils;

import paramfile.ParameterSet;
import paramfile.parameters.DoubleParameter;
import paramfile.parameters.IntegerParameter;
import paramfile.parameters.StringParameter;
import paramfile.parameters.ValidationFailedException;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import data.GrmReader;
import data.IndividualId;
import data.IndividualVariables;
import data.PcSubtractingGrmReader;

public class JacknifingLowMemMultipleHeRegressionEstimatorBin {
	private static final String GRMS = "grmlist";
	private static final String PHENOS = "phenos";
	private static final String COVARS = "covars";
	private static final String PREVALENCE = "prevalence";
	private static final String MODCOUNT = "modcount";
	private static final String SUBTRACTVECS = "subtractvecs";
	private static final String COHORTS = "crosscohorts";
	
	private static final int PARTITIONS = 200;
	private static final double NULLSENTINEL = -1*Double.MAX_VALUE;
	protected static final String RSCRIPTNAME = "getConvertedPheno.r";
	
	private static Logger logger = Logger.getLogger(JacknifingLowMemMultipleHeRegressionEstimatorBin.class.getName());
	static
	{
		logger.setLevel(Level.INFO);
	}
	private static ParameterSet getParameters(String filename) throws FileNotFoundException, IOException, ValidationFailedException {
		ParameterSet params = new ParameterSet();
		params.addRequiredParameter(new StringParameter(GRMS));
		params.addRequiredParameter(new StringParameter(PHENOS));
		params.addParameter(new StringParameter(SUBTRACTVECS));
		params.addParameter(new StringParameter(COVARS));
		params.addRequiredParameter(new DoubleParameter(PREVALENCE));
		params.addRequiredParameter(new IntegerParameter(MODCOUNT, Integer.MAX_VALUE));
		params.addParameter(new StringParameter(COHORTS));
		params.load(filename);
		return params;
	}
	
	private static String getTempFileName(String base)
	{
		String temp_base = "temp-"+base;
		String file_name = temp_base;
		int index = 1;
		while((new File(file_name)).exists())
		{
			file_name = temp_base + index;
			index++;
		}
		return file_name;
	}
	
	//starting with a random assignment of individuals to cohorts
	public static void main(String[] args) throws IOException, ClassNotFoundException, ValidationFailedException, InterruptedException, URISyntaxException
	{
		ParameterSet p = getParameters(args[0]);
		String grmlist = p.getStringValue(GRMS);
		String untransformed_phenofile = p.getStringValue(PHENOS);
		String covars = p.getStringValue(COVARS);
		if( covars == null )
		{
			covars = untransformed_phenofile;
		}
		double prevalence = p.getDoubleValue(PREVALENCE);
		String subtractvecs = p.getStringValue(SUBTRACTVECS);
		
		String multiplier_filename = getTempFileName("multiplier");
		String transformed_filename = getTempFileName("transformed");
		(new File(multiplier_filename)).deleteOnExit();
		(new File(transformed_filename)).deleteOnExit();
		
		String jar = JacknifingLowMemMultipleHeRegressionEstimatorBin.class.getProtectionDomain().getCodeSource().getLocation().toString();
		int dir_start = jar.indexOf('/');
		int dir_stop = jar.lastIndexOf('/');
		String dir = jar.substring(dir_start, dir_stop);
		String rscript = dir + "/" + RSCRIPTNAME;
		Process process = Runtime.getRuntime().exec("Rscript  " + rscript + " " + untransformed_phenofile + " " + prevalence + " " + covars + " " + transformed_filename + " " +  multiplier_filename);
		Scanner scanner = new Scanner(process.getInputStream());
		double total_var = -1;
		while(scanner.hasNextDouble())
		{
			total_var = scanner.nextDouble();
		}
		Scanner error_scanner = new Scanner(process.getErrorStream());
		while(error_scanner.hasNext())
		{
			System.err.println(error_scanner.nextLine());
		}
		int exit = process.waitFor();
		if( (exit != 0) || (total_var == -1) )
		{
			System.err.println("R transformation of phenotype appears to have failed with " + exit);
			
		}
		
		int modcount = p.getIntValue(MODCOUNT);
		logger.log(Level.INFO, "reading GRMS");
		List<GrmReader> l = getGrms(grmlist, subtractvecs);
		IndividualVariables ids = l.get(0).getIds();
		List<IndividualId> ordered_ids = ids.getOrderedIndividualIds();
		logger.log(Level.INFO, "reading phenotypes");
		IndividualVariables pheno_set = IndividualVariables.read(transformed_filename);
		double[] mult_array = null;
		if( multiplier_filename != null )
		{
			IndividualVariables multiplier_set = IndividualVariables.read(multiplier_filename);
			mult_array = multiplier_set.getArray(ordered_ids, 0, NULLSENTINEL);
		}
		List<String> cohort_names = new ArrayList<String>();
		logger.log(Level.INFO, "reading cohorts");
		int[] cohorts = getCohorts(ordered_ids, p.getStringValue(COHORTS), cohort_names);
		logger.log(Level.INFO, "producing partitions");
		int[] partitions = getPartitionAssignments(ordered_ids.size(), PARTITIONS);
		int[][] partition_mapping = new int[PARTITIONS][];
		for(int i = 0; i < PARTITIONS; i++)
		{
			partition_mapping[i] = new int[PARTITIONS];
		}
		int partition_pair_index = 0;
		for(int i = 0; i < PARTITIONS; i++)
		{
			for(int j = i; j < PARTITIONS; j++)
			{
				partition_mapping[i][j] = partition_pair_index;
				partition_mapping[j][i] = partition_pair_index;
				partition_pair_index++;
			}
		}
		
		double[] pheno_array = pheno_set.getArray(ordered_ids, 0, NULLSENTINEL);
		int partition_pairs = PARTITIONS*(PARTITIONS + 1)/2;
		DoubleMatrix2D[] dens = new DoubleMatrix2D[partition_pairs];
		DoubleMatrix2D[] nums = new DoubleMatrix2D[partition_pairs];
		for(int i = 0; i < dens.length; i++)
		{
			dens[i] = new DenseDoubleMatrix2D(l.size(), l.size());
			nums[i] = new DenseDoubleMatrix2D(l.size(), 1);
		}
		int count = 0;
		while(l.get(0).hasNext())
		{
			count += addNext(l, dens, nums, pheno_array, mult_array, partitions, partition_mapping, cohorts);
			if((count % modcount == 0))
			{
				System.out.print(count + " pairs analyzed\r");
			}
		}
		System.out.println();
		DoubleMatrix2D sum_den = sumAll(dens);
		DoubleMatrix2D sum_num = sumAll(nums);
//		print(sum_num);
//		print(sum_den);
//		System.out.println(total_var);
		double[] full_estimates = estimate(total_var, sum_num, sum_den);
		double[] ses = new double[full_estimates.length];
		double[][] jacknife_estimates = new double[full_estimates.length][];
		for(int i = 0; i < jacknife_estimates.length; i++)
		{
			jacknife_estimates[i] = new double[PARTITIONS - 1];
		}
		for(int index = 1; index < PARTITIONS; index++)
		{
			double[] result = subtractAndEstimate(total_var, sum_num, sum_den, nums, dens, index, partition_mapping, PARTITIONS);
			for(int component_index = 0; component_index < jacknife_estimates.length; component_index++)
			{
				jacknife_estimates[component_index][index - 1] = result[component_index];
			}
		}
		for(int i = 0; i < full_estimates.length; i++)
		{
			double variance = StatUtils.variance(jacknife_estimates[i])*ordered_ids.size();
			ses[i] = Math.sqrt(variance);
		}
		for(int i = 0; i < full_estimates.length - 1; i++)
		{
			System.out.println("V(G" + i + ")/Vp_L:\t" + full_estimates[i] + "\tS.E.:\t" + ses[i]);
		}
		System.out.println("Sum of V(G)/Vp_L:\t" + full_estimates[full_estimates.length - 1] + "\tS.E.:\t" + ses[full_estimates.length - 1]);
	}
	
	private static int[] getCohorts(List<IndividualId> ids, String cohorts_file, List<String> cohort_names) throws FileNotFoundException
	{
		int[] cohort_id = new int[ids.size()];
		if( cohorts_file == null )
		{
			for(int i = 0; i < cohort_id.length; i++)
			{
				cohort_names.add("" + i);
				cohort_id[i] = i;
			}
			return cohort_id;
		}
		Scanner scanner = new Scanner(new BufferedInputStream(new FileInputStream(cohorts_file)));
		Map<String, Integer> cohort_number = new HashMap<String, Integer>();
		int index = 0;
		while(scanner.hasNext())
		{
			String fid = scanner.next();
			String iid = scanner.next();
			IndividualId id = new IndividualId(fid, iid);
			String cohort = scanner.next();
			if( !cohort_number.containsKey(cohort) )
			{
				cohort_number.put(cohort, cohort_number.size());
				cohort_names.add(cohort);
			}
			if( !ids.get(index).equals(id) )
			{
				scanner.close();
				throw new RuntimeException(ids.get(index) + " is not the same as " + id);
			}
			cohort_id[index] = cohort_number.get(cohort);
			index++;
		}
		scanner.close();
		return cohort_id;
	}
	
	private static double[] subtractAndEstimate(double totalvar, DoubleMatrix2D sum_num, DoubleMatrix2D sum_den, DoubleMatrix2D[] nums, DoubleMatrix2D[] dens, int to_subtract, int[][] partition_mapping, int partitions) {
		DoubleMatrix2D copy_num = sum_num.copy();
		DoubleMatrix2D copy_den = sum_den.copy();
		for(int i = 0; i < partitions; i++)
		{
			int mat_index = partition_mapping[to_subtract][i];
			copy_num.assign(nums[mat_index], Functions.minus);
			copy_den.assign(dens[mat_index], Functions.minus);
		}
		return estimate(totalvar, copy_num, copy_den);
		
	}
	
	public static void print(DoubleMatrix2D mat)
	{
		for(int i = 0; i < mat.rows(); i++)
		{
			for(int j = 0; j < mat.columns(); j++)
			{
				System.out.print(mat.get(i, j) + "\t");
			}
			System.out.println();
		}
	}
	
	private static double[] estimate(double totalvar, DoubleMatrix2D sum_num, DoubleMatrix2D sum_den)
	{
		Algebra a = new Algebra();
		DoubleMatrix2D inv = a.inverse(sum_den);
		DoubleMatrix2D coeffs = a.mult(inv, sum_num);
		double[] results = new double[sum_num.rows() + 1];
		double sum = 0;
		for(int i = 0; i < coeffs.rows(); i++)
		{
			double cg = coeffs.getQuick(i, 0);
			cg /= (1 + totalvar);
			sum += cg;
			results[i] = cg;
		}
		results[results.length - 1] = sum;
		return results;
	}

	private static DoubleMatrix2D sumAll(DoubleMatrix2D [] to_sum)
	{
		DoubleMatrix2D sum = new DenseDoubleMatrix2D(to_sum[0].rows(), to_sum[0].columns());
		for(DoubleMatrix2D part : to_sum )
		{
			sum.assign(part, Functions.plus);
		}
		return sum;
	}
	
	private static int[] getPartitionAssignments(int n, int partitions) throws FileNotFoundException
	{
		//partition # of individuals get their own partition.
		//Everyone else is in the main partition
		Random rand = new Random();
		int[] assignment = new int[n];
		for(int i = 1; i < partitions; i++)
		{
			int individual = rand.nextInt(n);
			while( assignment[individual] != 0 )
			{
				individual = rand.nextInt(n);
			}
			assignment[individual] = i;
		}
//		int partition_size = assignment.length/partitions;
//		for(int i = 0; i < assignment.length; i++)
//		{
//			assignment[i] = (i % partition_size);
//		}
//		MathArrays.shuffle(assignment);
		return assignment;
	}
	
//	private static Map<Integer, Integer> produceFinal(List<IndividualId> ordered_grm_id, IndividualVariables phenotypes, IndividualVariables covariates)
//	{
//		Map<Integer, Integer> mapped = new HashMap<Integer, Integer>();
//		int index = -1;
//		int mapped_index = -1;
//		for(IndividualId id : ordered_grm_id)
//		{
//			index++;
//			if( !phenotypes.contains(id) )
//			{
//				continue;
//			}
//			if( !covariates.contains(id) )
//			{
//				continue;
//			}
//			mapped_index++;
//			mapped.put(index, mapped_index);
//		}
//		//all ids were kept, return null
//		if( mapped.size() == ordered_grm_id.size() )
//		{
//			return null;
//		}
//		return mapped;
//	}

	private static int addNext(List<GrmReader> l, DoubleMatrix2D[] dens, DoubleMatrix2D[] nums, double[] phenos, double[] multiplier, int[] parititons, int[][] partition_mapping, int[] cohorts) throws IOException {
		double[][] vals = new double[l.size()][];
		int index1 = 0, index2 = 0;
		boolean keep_looking = true;
		while(keep_looking)
		{
			if(!l.get(0).hasNext())
			{
				logger.log(Level.INFO, "no entries after " + l.get(0).getIndex1() + "," + l.get(0).getIndex2());
				return 0;
			}
			l.get(0).next();
			vals[0] = new double[1];
			vals[0][0] = l.get(0).getValue();
			index1 = l.get(0).getIndex1();
			index2 = l.get(0).getIndex2();
			for(int i = 1; i < l.size(); i++)
			{
				GrmReader s  = l.get(i);
				if(!s.hasNext())
				{
					throw new RuntimeException();
				}
				s.next();
				vals[i] = new double[1];
				vals[i][0] = s.getValue();
			}
			keep_looking = false;
			if(index1 == index2)
			{
				keep_looking = true;
			}
			else if(phenos[index1] == NULLSENTINEL) 
			{
				keep_looking = true;
				logger.log(Level.WARNING, "skipped pair because individual " + index1 + " has a null phenotype.");
			}
			else if(phenos[index2] == NULLSENTINEL) 
			{
				keep_looking = true;
				logger.log(Level.WARNING, "skipped pair because individual " + index2 + " has a null phenotype.");
			}
			else if((multiplier != null) && ((multiplier[index1] == NULLSENTINEL) || (multiplier[index2] == NULLSENTINEL)) )
			{
				keep_looking = true;
				logger.log(Level.WARNING, "skipped pair because individual " + index1 + " or " + index2 + " has a null multiplier.");
			}
			else if(cohorts[index1] == cohorts[index2])
			{
				keep_looking = true;
			}
		}
		if( multiplier != null )
		{
			for(int index = 0; index < vals.length; index++)
			{
				vals[index][0] *= multiplier[index1]*multiplier[index2];
			}
		}
		double yval = phenos[index1]*phenos[index2];
		int cohort_pair_index = partition_mapping[parititons[index1]][parititons[index2]];
		DoubleMatrix2D d2d = new DenseDoubleMatrix2D(vals);
		d2d.zMult(d2d, dens[cohort_pair_index], 1, 1, false, true);
		nums[cohort_pair_index].assign(d2d, Functions.plusMult(yval));
		return 1;
	}


	public static List<GrmReader> getGrms(String grmlist, String subtract_vecs) throws IOException
	{
		List<String> grm_prefixes = getStrings(grmlist);
		List<String> vec_prefixes = null;
		if( subtract_vecs != null )
		{
			vec_prefixes = getStrings(subtract_vecs);
			if( vec_prefixes.size() != grm_prefixes.size() )
			{
				throw new RuntimeException();
			}
		}
		List<GrmReader> l = new ArrayList<GrmReader>();
		IndividualVariables ids = null;
		for(int i = 0; i < grm_prefixes.size(); i++)
		{
			String prefix = grm_prefixes.get(i);
			logger.log(Level.INFO, "reading GRM with prefix " + prefix);
			String idfile = prefix + ".grm.id";
			String binfile = prefix + ".grm.bin";
			GrmReader reader = null;
			if( vec_prefixes != null )
			{
				String vec_prefix = vec_prefixes.get(i);
				String vecfile = vec_prefix + ".eigenvec";
				String valfile = vec_prefix + ".eigenval";
				logger.log(Level.INFO, "reading subtracting vectors with prefix " + vec_prefix);
				reader = new PcSubtractingGrmReader(binfile, idfile, vecfile, valfile);
			}
			else
			{
				reader = new GrmReader(binfile, idfile);
			}
			logger.log(Level.INFO, "read " + reader.getIds().getOrderedIndividualIds().size() + " ids from GRM");
			if( ids == null )
			{
				ids = reader.getIds();
			}
			else if( !reader.getIds().getOrderedIndividualIds().equals(ids.getOrderedIndividualIds()) )
			{
				logger.log(Level.SEVERE, "ids for grm " + prefix + " not identical to previous");
				throw new RuntimeException(idfile + " not identical to previous\n");
			}
			l.add(reader);
		}
		return l;
	}
	
	public static List<String> getStrings(String listfile) throws IOException
	{
		Scanner scanner = null;
		List<String> l;
		try {
			scanner = new Scanner(new BufferedInputStream(new FileInputStream(listfile)));
			
			l = new ArrayList<String>();
			while(scanner.hasNext())
			{
				String s = scanner.next();
				l.add(s);
			}
		} finally
		{
			if( scanner != null )
			{
				scanner.close();
			}
		}
		return l;
	}
	
	
}
