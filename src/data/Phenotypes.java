package data;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

public class Phenotypes {
	private static final String MISSING1 = "-9";
	private static final String MISSING2 = "NA";
	private Map<IndividualId, Double> m_phenotypes;
	private double m_mean_phenotype;
	
	
	private Phenotypes(Map<IndividualId, Double> phenotypes)
	{
		m_phenotypes = phenotypes;
		double sum = 0;
		int count = 0;
		for(Double val : phenotypes.values())
		{
			sum += val;
			count++;
		}
		m_mean_phenotype = sum/count;
	}
	
	public double getPhenotype(IndividualId id)
	{
		return m_phenotypes.get(id);
	}
	
	public double getMeanPhenotype()
	{
		return m_mean_phenotype;
	}
	
	public Set<IndividualId> getIndividualIds() {
		return m_phenotypes.keySet();
	}
	
	public Phenotypes adjust(IndividualVariables covars)
	{
		List<IndividualId> kept_ids = new ArrayList<IndividualId>(m_phenotypes.keySet());
		for(IndividualId id : covars.getIndividualIds())
		{
			if( !m_phenotypes.containsKey(id) )
			{
				continue;
			}
			boolean all_non_null = true;
			for(Double covariate_value : covars.getValues(id))
			{
				if( covariate_value == null )
				{
					all_non_null = false;
					break;
				}
			}
			if( !all_non_null )
			{
				continue;
			}
			kept_ids.add(id);
		}
		
		double[] y = new double[kept_ids.size()];
		double[][] x = new double[kept_ids.size()][];
		for(int i = 0; i < kept_ids.size(); i++)
		{
			y[i] = getPhenotype(kept_ids.get(i));
			x[i] = covars.getArray(kept_ids.get(i));
		}
		OLSMultipleLinearRegression reg = new OLSMultipleLinearRegression();
		reg.newSampleData(y, x);
		double[] resid = reg.estimateResiduals();
		Map<IndividualId, Double> phenos = new HashMap<IndividualId, Double>();
		for(int i = 0; i < resid.length; i++)
		{
			phenos.put(kept_ids.get(i), resid[i]);
		}
		return new Phenotypes(phenos);
	}
	
	public static Phenotypes read(String phenofile) throws FileNotFoundException
	{
		Map<IndividualId, Double> phenotypes = new HashMap<IndividualId, Double>();
		Scanner scanner = new Scanner(new BufferedInputStream(new FileInputStream(phenofile)));
		while(scanner.hasNext())
		{
			String fid = scanner.next();
			String iid = scanner.next();
			String string_phenotype = scanner.next();
			if( string_phenotype.equals(MISSING1) || string_phenotype.equals(MISSING2))
			{
				continue;
			}
			double phenotype = Double.parseDouble(string_phenotype);
			IndividualId id = new IndividualId(fid, iid);
			phenotypes.put(id, phenotype);
		}
		scanner.close();
		return new Phenotypes(phenotypes);
	}
	
}
