package image_processing;

import java.math.BigInteger;
import java.security.NoSuchAlgorithmException;
import java.sql.*;
import java.io.*;
import java.net.*;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Random;

public class DBUtil {

	private static InputStream is;
	private static OutputStream os;
	private static BufferedWriter bw;
	private static BufferedReader br;
	private static DataInputStream dr;
	private static DataOutputStream dw;
	private static final String URL = "jdbc:mysql://127.0.0.1:3306/image library?useUnicode=true&characterEncoding=UTF-8";
	private static final String USER = "root";
	private static final String PASSWORD = "20010609";

	private static Connection conn = null;

	private static ServerSocket ss;
	private static Socket Server;

	private static int Kappa=128,K=448,M=100000,N=10000,L=96;


	public static int[][] node_id=new int[15][];


	static {
		try {
			// 1.������������
			Class.forName("com.mysql.jdbc.Driver");
			// 2.������ݿ������
			conn = DriverManager.getConnection(URL, USER, PASSWORD);
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}

	public static void BuildServer() {
		try {
			//3.服务器端
			ss = new ServerSocket(8888);
			System.out.println("启动服务器....");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	public static void AcceptClient() {
		try {
			Server = ss.accept();
			System.out.println("客户端:"+Server.getInetAddress().getLocalHost()+"已连接到服务器");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void TestCommunication() throws IOException {
		try {
			is = Server.getInputStream();
			os = Server.getOutputStream();
			br = new BufferedReader(new InputStreamReader(is));
			//读取客户端发送来的消息
			String mess = br.readLine();
			System.out.println("客户端："+mess);
			bw = new BufferedWriter(new OutputStreamWriter(os));
			bw.write(mess+"\n");
			bw.flush();
		}catch (IOException e) {
			e.printStackTrace();
		}
		//br.close();
		//bw.close();
	}

	public static Connection getConnection() {
		return conn;
	}

	private static LinkedList<Integer> Selectionindexes=new LinkedList <Integer>();
	private static byte[][] uArray=new byte[K][];
	private static byte[][] vArray=new byte[K][];
	private static byte[][] djArray=new byte[M][];

	private static byte[][] diArray=new  byte[K][];

	private static byte[][] t0Array=new  byte[K][];

	private static byte[][] t1Array=new  byte[K][];

	private static byte[][] t0jArray=new  byte[M][];
	public static LinkedList<byte []> CodeWords;

	public static void Generate_OT() throws NoSuchAlgorithmException, IOException {
		CodeWords=CodeFunction.ZeroOne(Kappa);
		PreProcessing_Receiver(K,2,Kappa);
		PictureShapeSimilarityContrast.PreProcessing_Sender(K,2,Kappa,Kappa);
		BaseDDH_OT(K,2,Kappa);
		OpenDataIO();
		PictureShapeSimilarityContrast.OpenDataIO();
		OT_Extension(K,2,Kappa,Kappa);
	}
	public static void PreProcessing_Receiver(int m, int n, int k) throws NoSuchAlgorithmException {
		//random choice r=(r1,r_2,...)
		for(int j=0; j<m; j++) {
			Random rn = new Random();
			Integer randomNum = rn.nextInt(n);
			Selectionindexes.add(randomNum);
		}

		//Set matrix U and V
		for(int i=0; i < k; i++) {
			byte[] bytes = new byte[k/8];
			new Random().nextBytes(bytes);
			uArray[i] = bytes;
		}
		for(int i=0; i < k; i++) {
			byte[] bytes = new byte[k/8];
			new Random().nextBytes(bytes);
			vArray[i] = bytes;
		}

		//Set djArray: Code Matrix
		for (int j= 0; j<m; j++) {
			djArray[j]= CodeWords.get(Selectionindexes.get(j));
		}

		//Set diArray: Code Matrix Transposition
		int height = m;
		int width = k;
		int[][] array = new int[width][height];
		for (int j=0; j<width; j++) {
			for (int i = 0; i < height; i++) {
				array[j][i] = getBit(djArray[i], j);
			}
		}
		for (int i = 0; i < k; i++) {
			diArray[i] = intArrayToByteArray(array[i]);
		}

		//Set TArray:T0 and T1
		for (int i = 0; i < k; i++) {
			t1Array[i] = RandomOracles.G(m/8,vArray[i]);
			t0Array[i] = xorByteArrays(t1Array[i],diArray[i]);
		}

		//Set T0JArray: Transposition of T0
		width = m;
		height = k;
		array = new int[width][height];
		for (int i = 0; i <width; i++) {
			for (int j = 0; j < height; j++) {
				array[i][j] = getBit(t0Array[j], i);
			}
		}
		for (int i = 0; i < m; i++) {
			t0jArray[i] = intArrayToByteArray(array[i]);
		}
	}

	public static void BaseDDH_OT(int m,int n,int k) throws IOException {
		System.out.println();
		System.out.println("\t\tBaseDDH_OT phase begins now...");
		System.out.println();

		for (int i = 0; i < k; i++) {
			PictureShapeSimilarityContrast.BaseDDH_ithOT_GenerateKeys(i);

			BigInteger keys[]=new BigInteger[2];
			keys[0]=new BigInteger(br.readLine());
			keys[1]=new BigInteger(br.readLine());

			BigInteger x0,x1;
			x0 = UnsignedBigInteger.fromUnsignedByteArray(uArray[i]);
			x1 = UnsignedBigInteger.fromUnsignedByteArray(vArray[i]);
			BigInteger r=ElGamalEncryption.getRandom();
			BigInteger[] cipher0=ElGamalEncryption.encrypt(keys[0], r, x0);
			BigInteger[] cipher1=ElGamalEncryption.encrypt(keys[1], r, x1);


			bw.write(cipher0[0].toString());// g<sup>r</sup> mod p
			bw.newLine();
			bw.write(cipher0[1].toString());//encrypted message
			bw.newLine();
			bw.write(cipher1[0].toString());// g<sup>r</sup> mod p
			bw.newLine();
			bw.write(cipher1[1].toString());//encrypted message
			bw.newLine();
			bw.flush();

			PictureShapeSimilarityContrast.BaseDDH_ithOT_ReceiveMsg(i,k);
		}

		System.out.println();
		System.out.println("\t\tBaseDDH_OT phase ends now...");
		System.out.println();
	}

	private static byte[][] wArray=new byte[K][];
	private static LinkedList<byte[]> Y_Recieved=new LinkedList<byte[]>();
	private static LinkedList <byte[]> Z=new LinkedList<byte[]>();
	public static void OT_Extension(int m,int n,int k,int l) throws IOException, NoSuchAlgorithmException {
		//br.close();
		//bw.close();
		//Set wArray:t0 xor G(u)
		for (int i = 0; i < k; i++) {
			wArray[i]=(xorByteArrays(RandomOracles.G(m/8,uArray[i]),t0Array[i]));
		}

		{
			System.out.println();
			System.out.println("\t\tSending wArray in IKNP's Receiver begins now...");
			System.out.println();

			for (int i = 0; i < k; i++) {
				dw.writeInt(wArray[i].length);
				dw.write(wArray[i]);
				dw.flush();
				PictureShapeSimilarityContrast.W_TransferReceiver(i);
			}

			System.out.println();
			System.out.println("\t\tSending wArray IKNP's Receiver ends now...");
			System.out.println();
		}

		PictureShapeSimilarityContrast.setQ_QjArray_andY(m,n,k,l);

		{
			System.out.println();
			System.out.println("\t\t Recieving Y IKNP's Receiver begins now ...");
			System.out.println();

			for (int j = 0; j < m; j++) {

				for (int nn = 0; nn < n; nn++) {
					byte[] y = new byte[l/8];
					PictureShapeSimilarityContrast.Y_TransferSender(n,j,nn);
					dr.readFully(y,0,l/8);
					Y_Recieved.add(y);
				}
			}

			int rj;
			byte[] hash;
			byte[] yj;
			for (int j = 0; j < m; j++) {
				rj=Selectionindexes.get(j);
				yj=Y_Recieved.get(j*n+rj);
				hash=RandomOracles.G(l/8, t0jArray[j]);
				Z.add(xorByteArrays(hash,yj));
			}

			/*System.out.println(" \n Selectionindexes:" + Selectionindexes);
			System.out.println("\n final result in Reciever :" +Z);*/

			System.out.println();
			System.out.println("\t\t Recieving Y IKNP's Receiver ends now ...");
			System.out.println();
		}
	}


	private static byte[] sArray;
	public static void PreProcessing2_Sender(int m,int n,int k,int l) {
		//Set sArray
		int[] bits=new int[k];
		for(int i=0;i<k;i++)
			bits[i]=Selectionindexes.get(i);
		sArray=intArrayToByteArray(bits);
	}

	private static byte[][] wArray_Recieved=new byte[K][];
	public static void W_TransferReceiver(int i) throws IOException {
		int size = dr.readInt();
		byte [] w=new byte[size];
		dr.readFully(w,0,size);
		wArray_Recieved[i]=w;
	}

	private static byte[][] qArray=new  byte[K][];
	private static byte[][] qjArray=new  byte[M][];
	public static void setQ_QjArray(int m, int n, int k, int l) throws NoSuchAlgorithmException {
		System.out.println();
		System.out.println("\t\tSet Q matrix in OPRF's Sender begins now...");
		System.out.println();

		for (int i = 0; i < k; i++) {
			if (Selectionindexes.get(i) == 0)  {
				qArray[i] =xorByteArrays(RandomOracles.G(m/8,Z.get(i)),wArray_Recieved[i]);
			}
			else if (Selectionindexes.get(i) == 1) {
				qArray[i] = RandomOracles.G(m/8,Z.get(i));
			}
		}

		int width = m;
		int height = k;
		int[][] array = new int[width][height];
		// Make a two-dimensional int array to put each bit of the byte arrays
		// Transposition of the byte array to the int array
		for (int i = 0; i<width ; i++) {
			for (int j = 0; j < height; j++) {
				array[i][j] = getBit(qArray[j], i);
			}
		}
		// Convert each int array to byte array and insert it in t0jArray
		for (int i = 0; i < m; i++) {
			qjArray[i] = intArrayToByteArray(array[i]);
		}

		System.out.println();
		System.out.println("\t\tSet Q matrix in OPRF's Sender ends now...");
		System.out.println();
	}
	private static byte[] Y;
	public static void Set_Y(int i,int dis) throws SQLException, NoSuchAlgorithmException {
		Y=new byte[(node_id[i].length-1)*L/8];
		Connection con = getConnection();
		Statement stmt = con.createStatement();
		ResultSet rs = stmt.executeQuery("select * from shape_kd_tree");
		byte[]and;
		byte[]xor;
		for(int nn=0;nn<node_id[i].length-1;nn++) {
			rs.absolute(node_id[i][nn+1]);
			byte[] SPLIT_DIM,FV;
			int split_dim=rs.getInt("SPLIT_DIM");
			SPLIT_DIM=ByteUtils.intToByteArray(split_dim);
			if(split_dim==0) FV=ByteUtils.double2Bytes(rs.getDouble("FV0"));
			else if(split_dim==1) FV=ByteUtils.double2Bytes(rs.getDouble("FV1"));
			else FV=ByteUtils.double2Bytes(rs.getDouble("FV2"));
			byte[] res=new byte[12];
			ByteUtils.ByteAdd(res,SPLIT_DIM,0,0,4);
			ByteUtils.ByteAdd(res,FV,4,0,8);
			and=PictureShapeSimilarityContrast.andByteArrays(PictureShapeSimilarityContrast.CodeWords.get((nn-dis+N)%N),sArray);
			xor=xorByteArrays(qjArray[PictureShapeSimilarityContrast.OT_Round],and);
			byte[] H_val=RandomOracles.G(12,xor);
			byte[] y_nn=xorByteArrays(res,H_val);
			ByteUtils.ByteAdd(Y,y_nn,nn*12,0,L/8);
		}
	}

	public static void Y_TransferSender() throws IOException {
		dw.write(Y);
		dw.flush();
	}

	public static void CloseBufferIO() throws IOException {
		br.close();
		bw.close();
	}
	public static void OpenDataIO() {
		dr=new DataInputStream(new BufferedInputStream(is));
		dw=new DataOutputStream(new BufferedOutputStream(os));
	}
	public static void CloseDataIO() throws IOException {
		dr.close();
		dw.close();
	}

	private static int getBit(byte[] data, int pos) {
		int posByte = pos/8;
		int posBit = pos%8;
		byte valByte = data[posByte];
		int valInt = valByte>>(8-(posBit+1)) & 0x0001;
		return valInt;
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


	public static byte[] xorByteArrays(byte[] a, byte[] b) {

		byte[] result = new byte[Math.min(a.length, b.length)];

		if (!(a.length == b.length)) {
			System.out.println("Lengths NOT equal");
			System.exit(3);
		}

		for (int i = 0; i < result.length; i++) {
			result[i] = (byte) (((int) a[i]) ^ ((int) b[i]));
		}
		return result;
	}

}
