
package image_processing;

public class ByteUtils {

  /**
   * Return a new byte array containing a sub-portion of the source array
   * @param srcBegin
   * The beginning index (inclusive)
   * @return The new, populated byte array
   **/
  public static byte[] subbytes(byte[] source, int srcBegin) {
    return subbytes(source, srcBegin, source.length);
  }
  /**
   * Return a new byte array containing a sub-portion of the source array
   * 
   * @param srcBegin
   *          The beginning index (inclusive)
   * @param srcEnd
   *          The ending index (exclusive)
   * @return The new, populated byte array
   **/
  public static byte[] subbytes(byte[] source, int srcBegin, int srcEnd) {
    byte destination[];

    destination = new byte[srcEnd - srcBegin];
    getBytes(source, srcBegin, srcEnd, destination, 0);

    return destination;
  }

  public static byte[] byteMerger(byte[] byte_1, byte[] byte_2){
    byte[] byte_3 = new byte[byte_1.length+byte_2.length];
    System.arraycopy(byte_1, 0, byte_3, 0, byte_1.length);
    System.arraycopy(byte_2, 0, byte_3, byte_1.length, byte_2.length);
    return byte_3;
  }
  public static byte[] byteIntMerger(byte[] byte_1, int i){
    byte[] byte_2=intToByteArray(i);
    return byteMerger(byte_1,byte_2);
  }
  public static int byteArrayToInt(byte[] b,int from) {
    return   b[3+from] & 0xFF |
            (b[2+from] & 0xFF) << 8 |
            (b[1+from] & 0xFF) << 16 |
            (b[from] & 0xFF) << 24;
  }

  public static byte[] intToByteArray(int a) {
    return new byte[] {
            (byte) ((a >> 24) & 0xFF),
            (byte) ((a >> 16) & 0xFF),
            (byte) ((a >> 8) & 0xFF),
            (byte) (a & 0xFF)
    };
  }


  /**
   * Copies bytes from the source byte array to the destination array
   * 
   * @param source
   *          The source array
   * @param srcBegin
   *          Index of the first source byte to copy
   * @param srcEnd
   *          Index after the last source byte to copy
   * @param destination
   *          The destination array
   * @param dstBegin
   *          The starting offset in the destination array
   */
  public static void getBytes(byte[] source, int srcBegin, int srcEnd, byte[] destination,
      int dstBegin) {
    System.arraycopy(source, srcBegin, destination, dstBegin, srcEnd - srcBegin);
  }

  private static int getBit(byte[] data, int pos) {
    int posByte = pos/8;
    int posBit = pos%8;
    byte valByte = data[posByte];
    int valInt = valByte>>(8-(posBit+1)) & 0x0001;
    return valInt;
  }

  public static void ByteAdd(byte[] A,byte[] B,int l1,int l2,int len) {
    for(int i=0;i<len;i++) {
      A[l1+i]=B[l2+i];
    }
  }

  public static byte[] double2Bytes(double d) {
    long value = Double.doubleToRawLongBits(d);
    byte[] byteRet = new byte[8];
    for (int i = 0; i < 8; i++) {
      byteRet[i] = (byte) ((value >> 8 * i) & 0xff);
    }
    return byteRet;
  }


  public static double bytes2Double(byte[] arr,int from) {
    long value = 0;
    for (int i = 0; i < 8; i++) {
      value |= ((long) (arr[i+from] & 0xff)) << (8 * i);
    }
    return Double.longBitsToDouble(value);
  }
  
}
