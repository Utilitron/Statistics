import java.math.BigInteger;
import java.util.List;

/**
 * Statistics Math
 */
public class StatsMath {
	public static float inverse(float val){ return 1-val; }
	
    public static BigInteger factorial(int n) {
    	BigInteger fact = BigInteger.valueOf(1); // this  will be the result
        for (int i = 1; i <= n; i++) {
        	fact = fact.multiply(BigInteger.valueOf(i));
        }
        return fact;
    }
	
	public static float mean(float[] list){ 
		float sum = 0;
	    for (int i = 0; i < list.length; i++) {
	        sum += list[i];
	    }
	    
	    float mean = sum / list.length;
	    
	    System.out.println("Mean: " + mean);
	    return mean;
	}
	
	public float median(float[] list){ 
	    int middle = list.length/2;
	    float median = 0;
	    if (list.length%2 == 1) {
	    	median = list[middle];
	    } else {
	    	median = (float) ((list[middle-1] + list[middle]) / 2.0);
	    }
	    
	    System.out.println("Median: " + median);
	    return median;  
	}
	
	public float mode(float[] list){ 
		float mode = 0;
	    int maxCount = 0;

	    for (int i = 0; i < list.length; ++i) {
	        int count = 0;
	        for (int j = 0; j < list.length; ++j) {
	            if (list[j] == list[i]) ++count;
	        }
	        if (count > maxCount) {
	            maxCount = count;
	            mode = list[i];
	        }
	    }

	    System.out.println("Mode: " + mode);
	    return mode;
	}
	
	public static float variance(float[] list){ 
		float mean = mean(list);
		float temp = 0;
        for(float a :list)
            temp += (mean-a)*(mean-a);
        
        float variance =  temp/list.length;
            
    	System.out.println("Variance: " + variance);
    	return variance;
	}
	
	public static float stdDev(float[] list){ 
		float stdDev = (float) Math.sqrt(variance(list));
		
		System.out.println("StdDev: " + stdDev);
    	return stdDev;
	}
	
	public static float confidenceInterval(float[] list, float magic){
		float n = (float) Math.sqrt(list.length);
		float confidenceInterval = (float) magic * (stdDev(list)/n);
		
		System.out.println("CI: " + confidenceInterval);
    	return confidenceInterval;
	}	
	
	
	
	public static float mean(List<Float> list){ 
		float sum = 0;
	    for (int i = 0; i < list.size(); i++) {
	        sum += list.get(i);
	    }
	    
	    float mean = sum / list.size();
	    
	    System.out.println("Mean: " + mean);
	    return mean;
	}
	
	public float median(List<Float> list){ 
	    int middle = list.size()/2;
	    float median = 0;
	    if (list.size()%2 == 1) {
	    	median = list.get(middle);
	    } else {
	    	median = (float) ((list.get(middle-1) + list.get(middle)) / 2.0);
	    }
	    
	    System.out.println("Median: " + median);
	    return median;  
	}
	
	public float mode(List<Float> list){ 
		float mode = 0;
	    int maxCount = 0;

	    for (int i = 0; i < list.size(); ++i) {
	        int count = 0;
	        for (int j = 0; j < list.size(); ++j) {
	            if (list.get(j) == list.get(i)) ++count;
	        }
	        if (count > maxCount) {
	            maxCount = count;
	            mode = list.get(i);
	        }
	    }

	    System.out.println("Mode: " + mode);
	    return mode;
	}
	
	public static float variance(List<Float> list){ 
		float mean = mean(list);
		float temp = 0;
        for(float a :list)
            temp += (mean-a)*(mean-a);
        
        float variance =  temp/list.size();
            
    	System.out.println("Variance: " + variance);
    	return variance;
	}
	
	public static float stdDev(List<Float> list){ 
		float stdDev = (float) Math.sqrt(variance(list));
		
		System.out.println("StdDev: " + stdDev);
    	return stdDev;
	}
	
	public static float confidenceInterval(List<Float> list, float magic) {
		float n = (float) Math.sqrt(list.size());
		float confidenceInterval = (float) magic * (stdDev(list)/n);
		
		System.out.println("CI: " + confidenceInterval);
    	return confidenceInterval;
	}
	

	
	public static float variance(float p_mean){ 
        float variance = p_mean*inverse(p_mean);
            
    	System.out.println("Variance: " + variance);
    	return variance;
	}
	
	public static float stdDev(float p_mean){ 
		float stdDev = (float) Math.sqrt(variance(p_mean));
		
		System.out.println("StdDev: " + stdDev);
    	return stdDev;
	}
	
	public static float confidenceInterval(float p_mean, int num, float magic){
		float n = (float) Math.sqrt(num);
		float confidenceInterval = (float) magic * (stdDev(p_mean)/n);
		
		System.out.println("CI: " + confidenceInterval);
    	return confidenceInterval;
	}		
	
	
	
	public static float binomialDistrobution(int n, int k){
		BigInteger fac_N = factorial(n);
		BigInteger fac_N_K = factorial(n-k);
		BigInteger fac_K = factorial(k);
		
		System.out.println("n!/(n-k)! k! :: " + n + "!/(" + n + "-" + k + ")!" + k + "!" + ":: " + fac_N + "/" + fac_N_K + "*" + fac_K);
		
		BigInteger b = fac_N.divide((fac_N_K.multiply(fac_K)));
		System.out.println("BinomialDistrobution: " + b.floatValue());

		return b.floatValue();
	}
	
	
	
	public float probability(float var, int num){
		float val = var;
		for (int i = 0; i < num; i++){
			val *= var;
		}
		
		return val; 
	}

	public float binomialProbability(int n, int k){
		return binomialProbability(n, k, 0.5f);
	}
	public static float binomialProbability(int n, int k, float probability){
		float bd = binomialDistrobution(n, k);
		float p_k = (float) Math.pow(probability,k);
		float p_nk = (float) Math.pow(inverse(probability),(n-k));
		
		System.out.println("(n!/(n-k)! k!) * (p^k) * ((1-p)^(n-k)) :: (" + n + "!/(" + n + "-" + k + ")!" + k + "!)" + "*(" + probability + "^" + k + ")*(" + inverse(probability) + "^(" + n + "-" + k +")) :: " + bd + "*" + p_k + "*" + p_nk);
		
		float p = bd * p_k * p_nk;
		System.out.println("Probability: " + p);
		return p;
	}
	
	/**
	 * Calculate the probability of a positive result given that
	 * p0=P(C)
	 * p1=P(Positive|C)
	 * p2=P(Negative|Not C)
	 * 
	 * @param probability
	 * @param sensitivity
	 * @param specificity
	 * @return positive
	 */
	public float positive(float probability, float sensitivity, float specificity){
		return probability*sensitivity+((inverse(probability))*(inverse(specificity)));
	}
	
	/**
	 * Calculate the probability of a negative result given that
	 * p0=P(C)
	 * p1=P(Positive|C)
	 * p2=P(Negative|Not C)
	 * 
	 * @param probability
	 * @param sensitivity
	 * @param specificity
	 * @return positive
	 */
	public float negative(float probability, float sensitivity, float specificity){
		return probability*inverse(sensitivity)/(probability*inverse(sensitivity)+inverse(probability) * specificity);
	}
	
	/**
	 * Return the probability of A conditioned on B given that 
	 * P(A)=p0, P(B|A)=p1, and P(Not B|Not A)=p2 
	 * 
	 * @param probability
	 * @param sensitivity
	 * @param specificity
	 * @return conditioned
	 */
	public float conditioned(float probability, float sensitivity, float specificity){
		return probability*sensitivity/positive(probability, sensitivity, specificity);
	}

	
	public static int[][] buildTruthTable(int n) {
	    int[][] rows = new int[(int) Math.pow(2,n)][n];

	    for (int i=0; i < rows.length; i++) {
	        for (int j= 0; j < rows[i].length; j++) {
	            rows[i][j] = (i/(int) Math.pow(2, j))%2;
	        }
	    }
	    return rows;
	}
	
	public static int countNumRows(int[][] truth){
		return truth.length;
	}
	
	public static int countRowsWithExact(int[][] truth, int val){
		int count = 0;
		
	    for (int i=0; i<truth.length; i++) {
	    	int rowVal = 0;
	        for (int j=0; j<truth[i].length; j++) {
	        	rowVal += truth[i][j];
	        }
	        
	        if (rowVal == val)
	        	count++;
	    }
		
		return count;
	}
	
	public static int countRowsWithAtLeast(int[][] truth, int val){
		int count = 0;
		
	    for (int i=0; i<truth.length; i++) {
	    	int rowVal = 0;
	        for (int j=0; j<truth[i].length; j++) {
	        	rowVal += truth[i][j];
	        }
	        
	        if (rowVal >= val)
	        	count++;
	    }
		
		return count;
	}
	
	public static void printTruthTable(int[][] truth){
	    for (int i=0; i<truth.length; i++) {
	        for (int j=0; j<truth[i].length; j++) {
	        	System.out.print(truth[i][j] + " ");
	        }
	        System.out.println("");
	    }
	}
}
