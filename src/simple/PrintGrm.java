package simple;

import java.io.IOException;

import data.GrmReader;

public class PrintGrm {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		GrmReader reader = new GrmReader(args[0] + ".grm.bin", args[0] + ".grm.id");
		while(reader.hasNext())
		{
			reader.next();
			System.out.println(reader.getIndex1() + "\t" + reader.getIndex2() + "\t" + reader.getValue());
		}
	}

}
