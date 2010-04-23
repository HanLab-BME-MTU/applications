// Here is an example of a concrete package (i.e. that implements the Package
// abstract class)
public class BioSensorsPackage extends Package {
	// Define the constructor
	public BioSensorsPackage(movieData owner) {
		super(owner, 'BioSensorsPackage')
	}
	
	// Following are the 3 abstract methods every package needs to implement.
	
	// Sanity check (See Package.m for more details)
	public void sanityCheck(boolean full) throws Exception {
		
		// If anything relative to the package itself needs to be checked before
		// we check the processes, do it here (throw an exception if needed)
		// ...
		
		// check the status of each process, the possible changed of their
		// parameters, and their dependencies with each other.
		checkProcesses_(full);
		
		// If anything relative to the package itself needs to be checked after
		// we've checked the processes, do it here (throw an exception if needed)
		// ...		
	}
	
	// This method returns the list of required processes. Each process is built
	// using the default contructor (See Package.m for more details)
	protected Process[] getDefaultProcessList_() {
		
		Process[] list = new Process[3];
		
		list[0] = new Name1Process(owner_);
		list[1] = new Name2Process(owner_);
		list[2] = new Name3Process(owner_);
		
		return list;
	}
	
	// This method return the dependency matrix of processes. It is a square
	// matrix whose size equals the number of required processes. (See Package.m
	// for more details)
	protected Matrix getProcessDependenyMatrix_() {
		return zeros(3);
	}
}