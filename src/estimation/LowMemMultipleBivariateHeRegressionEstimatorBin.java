package estimation;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import java.util.Scanner;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import data.GrmReader;
import data.IndividualId;
import data.IndividualVariables;

public class LowMemMultipleBivariateHeRegressionEstimatorBin {
	private static final String GRMS = "grmlist";
	private static final String PHENOS = "phenos";
	private static final String PHENOTYPE1 = "phenotype-index1";
	private static final String PHENOTYPE2 = "phenotype-index2";
	private static final String MULTIPLIER = "multiplier";
	private static final String TOTALVAR = "totalvar";
	private static final String MODCOUNT = "modcount";
	private static final String OUTFILE = "outfile";
	private static final String STOPCOUNT = "stopcount";
	private static final String COVARS = "covars";
	
	private static final double NULLSENTINEL = -1*Double.MAX_VALUE;
	
	public static void main(String[] args) throws IOException
	{
		String paramfile = args[0];
		Properties p = new Properties();
		p.load(new BufferedInputStream(new FileInputStream(paramfile)));
		String grmlist = p.getProperty(GRMS);
		String phenofile = p.getProperty(PHENOS);
		int phenotype_index1 = Integer.parseInt(p.getProperty(PHENOTYPE1));
		int phenotype_index2 = Integer.parseInt(p.getProperty(PHENOTYPE2));
		double total_var = Double.parseDouble(p.getProperty(TOTALVAR));
		
		int modcount = Integer.MAX_VALUE;
		if( p.getProperty(MODCOUNT) != null )
		{
			modcount = Integer.parseInt(p.getProperty(MODCOUNT));
		}
		int stopcount = -1; 
		if( p.getProperty(STOPCOUNT) != null )
		{
			stopcount = Integer.parseInt(p.getProperty(STOPCOUNT));
		}
		List<GrmReader> l = getGrms(grmlist);
		IndividualVariables ids = l.get(0).getIds();
		List<IndividualId> ordered_ids = ids.getOrderedIndividualIds();
		IndividualVariables pheno_set = IndividualVariables.read(phenofile);
		double[] mult_array = null;
		if( p.getProperty(MULTIPLIER) != null )
		{
			IndividualVariables multiplier_set = IndividualVariables.read(p.getProperty(MULTIPLIER));
			mult_array = multiplier_set.getArray(ordered_ids, 0, NULLSENTINEL);
		}
		if( p.getProperty(COVARS) != null )
		{
			IndividualVariables covariates = IndividualVariables.read(p.getProperty(COVARS));
			pheno_set = pheno_set.adjustAll(pheno_set.getOrderedIndividualIds(), covariates);
		}
		
		double[] pheno1_array = pheno_set.getArray(ordered_ids, phenotype_index1, NULLSENTINEL);
		double[] pheno2_array = pheno_set.getArray(ordered_ids, phenotype_index2, NULLSENTINEL);
		DoubleMatrix2D den = new DenseDoubleMatrix2D(l.size(), l.size() );
		DoubleMatrix2D num = new DenseDoubleMatrix2D(l.size(), 1);
		int count = 0;
		while(l.get(0).hasNext())
		{
			count += addNext(l, den, num, pheno1_array, pheno2_array, mult_array);
			if((count % modcount == 0))
			{
				System.err.println(count);
			}
			if( (stopcount > 0) && (count >= stopcount) )
			{
				break;
			}
		}
		if( p.getProperty(OUTFILE) != null )
		{
			String outfile = p.getProperty(OUTFILE);
			ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(outfile));
			oos.writeInt(count);
			oos.writeObject(den);
			oos.writeObject(num);
			oos.close();
		}
//		den.assign(Functions.div(count));
//		num.assign(Functions.div(count));
		
		Algebra a = new Algebra();
//		System.err.println("den is\n");
//		for(int i = 0; i < den.rows(); i++)
//		{
//			for(int j = 0; j < den.columns(); j++)
//			{
//				System.err.print(den.get(i, j) + "\t");
//			}
//			System.err.println();
//		}
//		System.err.println("num is\n");
//		for(int i = 0; i < num.rows(); i++)
//		{
//			for(int j = 0; j < num.columns(); j++)
//			{
//				System.out.print(num.get(i, j) + "\t");
//			}
//			System.out.println();
//		}
		DoubleMatrix2D inv = a.inverse(den);
		DoubleMatrix2D coeffs = a.mult(inv, num);
		for(int i = 0; i < coeffs.rows(); i++)
		{
			double cg = coeffs.getQuick(i, 0);
			cg /= (1 + total_var);
			System.out.println("C(G" + i + "):\t" + cg);
		}
	}
	
//	private static Map<Integer, Integer> produceFinal(List<IndividualId> ordered_grm_id, IndividualVariables phenotypes, IndividualVariables covariates, int pheno_index1, int pheno_index2)
//	{
//		Map<Integer, Integer> mapped = new HashMap<Integer, Integer>();
//		int index = -1;
//		int mapped_index = -1;
//		for(IndividualId id : ordered_grm_id)
//		{
//			index++;
//			if( (phenotypes.getValue(id, pheno_index1) == null) && (phenotypes.getValue(id, pheno_index2) == null)  )
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

	private static int addNext(List<GrmReader> l, DoubleMatrix2D den, DoubleMatrix2D num, double[] phenos1, double[] phenos2, double[] multiplier) throws IOException {
		double[][] vals = new double[l.size()][];
		int index1 = 0, index2 = 0;
		boolean keep_looking = true;
		double pheno1 = NULLSENTINEL, pheno2 = NULLSENTINEL;
		while(keep_looking)
		{
			if(!l.get(0).hasNext())
			{
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
				continue;
			}
			if( (phenos1[index1] == NULLSENTINEL) || (phenos2[index2] == NULLSENTINEL) )
			{
				if( (phenos2[index1] == NULLSENTINEL) || (phenos1[index2] == NULLSENTINEL) )
				{
					keep_looking = true;
					continue;
				}
				else
				{
					pheno1 = phenos2[index1];
					pheno2 = phenos1[index2];
				}
			}
			else
			{
				pheno1 = phenos1[index1];
				pheno2 = phenos2[index2];
			}
			if( (multiplier != null) && ((multiplier[index1] == NULLSENTINEL) || (multiplier[index2] == NULLSENTINEL)) )
			{
				keep_looking = true;
				continue;
			}
		}
		
		
		if( multiplier != null )
		{
			for(int index = 0; index < vals.length; index++)
			{
				vals[index][0] *= multiplier[index1]*multiplier[index2];
			}
		}
		if( (pheno1 == NULLSENTINEL) || (pheno2 == NULLSENTINEL) )
		{
			throw new RuntimeException();
		}
		double yval = pheno1*pheno2;
//		System.err.println(index1 + "\t" + index2 + "\t" + mapped_index1 + "\t" + mapped_index2 + "\t" + phenos1[mapped_index1] + "\t" +phenos2[mapped_index2] + "\t" + yval);
		DoubleMatrix2D d2d = new DenseDoubleMatrix2D(vals);
		d2d.zMult(d2d, den, 1, 1, false, true);
		num.assign(d2d, Functions.plusMult(yval));
		return 1;
	}


	public static List<GrmReader> getGrms(String grmlist) throws IOException
	{
		Scanner scanner = null;
		IndividualVariables ids = null;
		List<GrmReader> l;
		try {
			scanner = new Scanner(new BufferedInputStream(new FileInputStream(grmlist)));
			l = new ArrayList<GrmReader>();
			while(scanner.hasNext())
			{
				String prefix = scanner.next();
				String idfile = prefix + ".grm.id";
				String binfile = prefix + ".grm.bin";
				GrmReader reader = new GrmReader(binfile, idfile);
				if( ids == null )
				{
					ids = reader.getIds();
				}
				else if( !reader.getIds().getOrderedIndividualIds().equals(ids.getOrderedIndividualIds()) )
				{
					throw new RuntimeException(idfile + " not identical to previous\n");
				}
				l.add(reader);
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
