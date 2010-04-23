// Here is an example of a concrete package (i.e. that implements the Package
// abstract class)
public class BioSensorsPackage extends Package {
	// Define the constructor
	public BioSensorsPackage(movieData owner) {
		super(owner, 'BioSensorsPackage')
	}
	
	// Following are the 3 abstract methods every package needs to implement.
	
	// Sanity check
	public void sanityCheck(boolean full) throws Exception {
		// check the status of each process and their dependencies with
		// each other.
		checkProcessDependencies_(full);
		
		// If needed, we can add specific checks below:
		
		// if (full) {
		// ...
		// } else {
		// ...
		// }
	}
	
	// This method returns the list of required processes. Each process is built
	// using the default contructor
	protected Process[] getDefaultProcessList_() {
		
		Process[] list = new Process[3];
		
		list[0] = new Name1Process(owner_);
		list[1] = new Name2Process(owner_);
		list[2] = new Name3Process(owner_);
		
		return list;
	}
	
	// This method return the dependency matrix of processes. It is a square
	// matrix whose size equals the number of required processes.
	protected Matrix getProcessDependenyMatrix_() {
		return eye(3);
	}
}