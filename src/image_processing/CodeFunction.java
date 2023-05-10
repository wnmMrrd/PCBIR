package image_processing;

import java.util.LinkedList;
import javax.crypto.Cipher;
import javax.crypto.spec.IvParameterSpec;
import javax.crypto.spec.SecretKeySpec;
import org.apache.commons.codec.binary.Base64;
public class CodeFunction {
    private static final String key="hj7x89H$yuBI0456";
    private static final String iv ="NIfb&95GUY86Gfgh";

    private static int K=448;

    public static String encryptAES(String data) throws Exception {
        try {
            Cipher cipher = Cipher.getInstance("AES/CBC/NoPadding");
            int blockSize = cipher.getBlockSize();
            byte[] dataBytes = data.getBytes();
            int plaintextLength = dataBytes.length;

            if (plaintextLength % blockSize != 0) {
                plaintextLength = plaintextLength + (blockSize - (plaintextLength % blockSize));
            }

            byte[] plaintext = new byte[plaintextLength];
            System.arraycopy(dataBytes, 0, plaintext, 0, dataBytes.length);

            SecretKeySpec keyspec = new SecretKeySpec(key.getBytes(), "AES");
            IvParameterSpec ivspec = new IvParameterSpec(iv.getBytes());  // CBC模式，需要一个向量iv，可增加加密算法的强度

            cipher.init(Cipher.ENCRYPT_MODE, keyspec, ivspec);
            byte[] encrypted = cipher.doFinal(plaintext);

            return encode(encrypted).trim(); // BASE64做转码。

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }
    public static String encode(byte[] byteArray) {
        return new String(new Base64().encode(byteArray));
    }

    public static LinkedList<byte []> PRC(int n) throws Exception {
        LinkedList<byte []> C= new LinkedList<byte []>();
        for(int i=0;i<n;i++) {
            String i_str=i+"";
            byte[] res=new byte[K/8];
            byte[] res1 =encryptAES("1"+i_str).getBytes();
            ByteUtils.ByteAdd(res,res1,0,0,16);
            byte[] res2 =encryptAES("2"+i_str).getBytes();
            ByteUtils.ByteAdd(res,res2,16,0,16);
            byte[] res3 =encryptAES("3"+i_str).getBytes();
            ByteUtils.ByteAdd(res,res3,32,0,16);
            byte[] res4 =encryptAES("4"+i_str).getBytes();
            ByteUtils.ByteAdd(res,res4,48,0,8);
            C.add(res);
        }
        return C;
    }
    public static LinkedList<byte []> Walsh_Hadamard(int k) {

        LinkedList<byte []> C= new LinkedList<byte []>();
        int N = k;
        int[] res=new int[N];
        boolean[][] H = new boolean[N][N];
        H[0][0] = true;
        for(int n = 1; n < N; n += n)
        {
            for(int i = 0; i < n; i++)
                for(int j = 0; j < n; j++)
                {
                    H[i+n][j] = H[i][j];
                    H[i][j+n] = H[i][j];
                    H[i+n][j+n] = !H[i][j];
                }
        }

        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                if(H[i][j]) res[j]=1;

                else res[j]=0;

            }

            C.add(intArrayToByteArray(res));

        }

        return C;

    }

    public static LinkedList<byte []> ZeroOne(int k) {
        LinkedList<byte []> C= new LinkedList<byte []>();
        int N = k;
        int[] res=new int[N];

        for(int i = 0; i < 2; i++)
        {
            for(int j = 0; j < N; j++)
            {
                if(i>0) res[j]=1;

                else res[j]=0;

            }

            C.add(intArrayToByteArray(res));

        }

        return C;
    }


    public static byte[] intArrayToByteArray(int[] bits) {

        // If the condition isn't satisfied, an AssertionError will be thrown.
        // The length MUST be divisible by 8.
        assert bits.length % 8 == 0;
        byte[] bytes = new byte[bits.length / 8];

        for (int i = 0; i < bytes.length; i++) {
            int b = 0;
            for (int j = 0; j < 8; j++)
                b = (b << 1) + bits[i * 8 + j];
            bytes[i] = (byte) b;
        }
        return bytes;
    }
}
