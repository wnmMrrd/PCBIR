
package image_processing;

import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

public class RandomOracles {

	// Method to create  Random Oracle G
	
    public static byte[] G (int offset,byte[] key) throws NoSuchAlgorithmException{
		byte[] res=new byte[offset];
		for(int i=0;i<offset/64;i++) {
			MessageDigest md = MessageDigest.getInstance("SHA-512");
			md.update(ByteUtils.byteIntMerger(key,i));
			byte[] encrypted = md.digest();
			ByteUtils.ByteAdd(res,encrypted,i*64,0,64);
		}
		if(offset%64!=0) {
			MessageDigest md = MessageDigest.getInstance("SHA-512");
			md.update(ByteUtils.byteIntMerger(key,offset/64));
			byte[] encrypted = md.digest();
			ByteUtils.ByteAdd(res,encrypted,(offset/64)*64,0,offset%64);
		}
		return res;
	}
    
    // Method to create Random Oracle H 
    
    public static int H (byte[] convertMe) throws NoSuchAlgorithmException {
    	
    	MessageDigest md = MessageDigest.getInstance("SHA-512");
    	md.update(convertMe); 
		byte[] encrypted = md.digest();
        return getBit(encrypted,0);
    }
    
    // Method to read the value of a bit (0 or 1) of a byte array
    
    private static int getBit(byte[] data, int pos) {
	      int posByte = pos/8; 
	      int posBit = pos%8;
	      byte valByte = data[posByte];
	      int valInt = valByte>>(8-(posBit+1)) & 0x0001;
	      return valInt;
	   }
}