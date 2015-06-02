package data;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Scanner;
import java.util.Set;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

public class IndividualVariables {
	private static final String MISSING1 = "-9";
	private static final String MISSING2 = "NA";
	private LinkedHashMap<IndividualId, List<Double>> m_values;
	private int m_variable_count;
	
	
	private IndividualVariables(LinkedHashMap<IndividualId, List<Double>> covars, int variable_count)
	{
		m_values = covars;
		m_variable_count = variable_count;
	}
	
	public List<Double> getValues(IndividualId id)
	{
		return m_values.get(id);
	}
	
	public Set<IndividualId> getIndividualIds() {
		return m_values.keySet();
	}
	public List<IndividualId> getOrderedIndividualIds() {
		return new ArrayList<IndividualId>(m_values.keySet());
	}
	
	
	public static IndividualVariables read(String file) throws FileNotFoundException
	{
		LinkedHashMap<IndividualId, List<Double>> covars = new LinkedHashMap<IndividualId, List<Double>>();
		int variable_count;
		Scanner scanner = null;
		try {
			scanner = new Scanner(new BufferedInputStream(new FileInputStream(file)));
			variable_count = -1;
			while(scanner.hasNext())
			{
				String line = scanner.nextLine();
				Scanner line_scanner = new Scanner(line);
				String fid = line_scanner.next();
				String iid = line_scanner.next();
				List<Double> values = new ArrayList<Double>();
				while(line_scanner.hasNext())
				{
					String string_covar = line_scanner.next();
					if( string_covar.equals(MISSING1) || string_covar.equals(MISSING2))
					{
						values.add(null);
					}
					else
					{
						double covar = Double.parseDouble(string_covar);
						values.add(covar);
					}
				}
				line_scanner.close();
				IndividualId id = new IndividualId(fid, iid);
				if( variable_count == -1 )
				{
					variable_count = values.size();
				}
				else if( values.size() != variable_count )
				{
					throw new RuntimeException("number of fields for " + id + "  is " + values.size() + " inconsistent with the prior count of " + variable_count);
				}
				covars.put(id, values);
			}
		} finally
		{
			if( scanner != null )
			{
				scanner.close();
			}
		}
		return new IndividualVariables(covars, variable_count);
	}

	public double[] getArray(IndividualId id) {
		List<Double> values = getValues(id);
		double[] array = new double[values.size()];
		for(int i = 0; i < array.length; i++)
		{
			array[i] = values.get(i);
		}
		return array;
	}
	
	public double[] getArray(int variable_index) {
		double[] array = new double[m_values.size()];
		int i = 0;
		for(IndividualId id : m_values.keySet())
		{
			array[i] = m_values.get(id).get(variable_index);
			i++;
		}
		return array;
	}
	
	public double[] getArray(List<IndividualId> ids, int variable_index, double null_sentinel) {
		double[] array = new double[ids.size()];
		int i = 0;
		for(IndividualId id : ids)
		{
			Double v = null;
			if( contains(id) )
			{
				v = m_values.get(id).get(variable_index);
			}
			if( v == null )
			{
				array[i] = null_sentinel;
			}
			else
			{
				array[i] = v;
			}
			i++;
		}
		return array;
	}
	
	public double[] getArray(List<IndividualId> ids, int variable_index) {
		double[] array = new double[ids.size()];
		int i = 0;
		for(IndividualId id : ids)
		{
			Double v = null;
			if( contains(id) )
			{
				v = m_values.get(id).get(variable_index);
			}
			if( v == null )
			{
				throw new NullPointerException();
			}
			else
			{
				array[i] = v;
			}
			i++;
		}
		return array;
	}
	
	public Double getValue(IndividualId id, int variable_index) {
		List<Double> values = getValues(id);
		return values.get(variable_index);
	}
	private void setValue(IndividualId id, int variable_index, Double value) {
		List<Double> values = getValues(id);
		values.set(variable_index, value);
	}
	
	public IndividualVariables adjustAll(List<IndividualId> to_keep, IndividualVariables covars)
	{
		LinkedHashMap<IndividualId, List<Double>> values = new LinkedHashMap<IndividualId, List<Double>>();
		for(IndividualId id : to_keep)
		{
			values.put(id, new ArrayList<Double>());
		}
		
		for(int current_variable = 0; current_variable < m_variable_count; current_variable++)
		{
			Set<IndividualId> kept = new LinkedHashSet<IndividualId>();
			Set<IndividualId> removed = new LinkedHashSet<IndividualId>();
			for(IndividualId id : values.keySet())
			{
				boolean all_non_null = true;
				if( m_values.get(id).get(current_variable) == null )
				{
					all_non_null = false;
				}
				else
				{
					for(Double covariate_value : covars.getValues(id))
					{
						if( covariate_value == null )
						{
							all_non_null = false;
							break;
						}
					}
				}
				if( !all_non_null )
				{
					removed.add(id);
				}
				else
				{
					kept.add(id);
				}
			}
			
			double[] y = new double[kept.size()];
			double[][] x = new double[kept.size()][];
			int array_index = 0;
			for(IndividualId id : kept)
			{
				y[array_index] = getValues(id).get(current_variable);
				x[array_index] = covars.getArray(id);
				array_index++;
			}
			OLSMultipleLinearRegression reg = new OLSMultipleLinearRegression();
			reg.newSampleData(y, x);
			double[] resid = reg.estimateResiduals();
			array_index = 0;
			for(IndividualId id : kept)
			{
				values.get(id).add(resid[array_index]);
				array_index++;
			}
			for(IndividualId id : removed)
			{
				values.get(id).add(null);
			}
		}
		
		return new IndividualVariables(values, m_variable_count);
	}

	public int getVariableCount() {
		return m_variable_count;
	}

	public void standardizeAll() {
		for(int i = 0; i < getVariableCount(); i++)
		{
			standardize(i);
		}
	}
	
	public void standardize(int variable_index)
	{
		double sum = 0;
		double sum_sq = 0;
		int count = 0;
		for(IndividualId id : getIndividualIds())
		{
			Double value = getValue(id, variable_index);
			if( value == null )
			{
				continue;
			}
			sum += value;
			sum_sq += value*value;
			count++;
		}
		double mean = sum/count;
		double sd = Math.sqrt(sum_sq/count - mean*mean);
		for(IndividualId id : getIndividualIds())
		{
			Double value = getValue(id, variable_index);
			if( value == null )
			{
				continue;
			}
			value = (value - mean)/sd;
			setValue(id, variable_index, value);
		}
	}

	public boolean contains(IndividualId id) {
		return m_values.containsKey(id);
	}
}
