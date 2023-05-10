package image_processing;

import java.math.BigInteger;
import java.security.NoSuchAlgorithmException;
import java.util.LinkedList;
import java.util.Random;
import javafx.util.Pair;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Arrays;
import java.io.*;
import java.net.*;

public class PictureShapeSimilarityContrast {

	private static InputStream is;
	private static OutputStream os;
	private static BufferedWriter bw;
	private static BufferedReader br;
	private static DataInputStream dr;
	private static DataOutputStream dw;
	private static double[] OriginalImageShape = new double[3];
	private static double[] ComparisonImageShape = new double[3];
	private static float ResultSet[];
	private static int ResultIDSet[];
	private static String ResultPathSet[];

	private static int node_cnt,Ls[],Rs[],Fa[],Ls_index[],level_size[],split_dim[];

	private static int minimal_id,Tm;

	private static double Fv[][];

	//private

	private static double minimal_distance,Fv_u[],ans_Fv[];

	private static double Original_fv[][],tmp_fv[][];

	private static Socket Client;

	private static int Kappa=128,K=448,M=100000,N=10000,L=96;

	/*
	 * 测试代码所用
	 */
	private static int mx_depth=0;
	public static void main(String[] args) throws Exception {

		try{
			DBUtil.BuildServer();
			Client = new Socket("127.0.0.1",8888);
			DBUtil.AcceptClient();

			//构建IO
			is = Client.getInputStream();
			os = Client.getOutputStream();

			bw = new BufferedWriter(new OutputStreamWriter(os));
			//向服务器端发送一条消息
			bw.write("测试客户端和服务器通信，服务器接收到消息返回到客户端\n");
			bw.flush();
			DBUtil.TestCommunication();

			//读取服务器返回的消息
			br = new BufferedReader(new InputStreamReader(is));
			String mess = br.readLine();
			System.out.println("服务器："+mess);
		} catch (UnknownHostException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		//CloseBufferIO();

		Generate_OT();

		int count=new GetRecordCount().getCount();
		prepare_testdata(count);
		final_test(count);
		//System.out.println(mx_depth);
		/*CloseBufferIO();
		CloseDataIO();
		DBUtil.CloseBufferIO();
		DBUtil.OpenDataIO();*/
	}
	// /Users/yangheng/Desktop/CBIR-master/datafiles/101_ObjectCategories/101_ObjectCategories/snoopy/image_0001.jpg

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
	private static int[] dep_to_node;
	public static void prepare_testdata(int count) throws SQLException {
		node_cnt=0;
		Ls=new int[(int)1e6];
		Rs=new int[(int)1e6];
		Fa=new int[(int)1e6];
		dep_to_node=new int[(int)1e6];
		Ls_index=new int[(int)1e6];
		level_size=new int[(int)1e6];
		split_dim=new int[(int)1e6];
		Fv=new double[3][(int)1e6];

		Original_fv=new double[3][count+1];
		tmp_fv=new double[3][count+1];

		Connection con = DBUtil.getConnection();
		Statement stmt = con.createStatement();
		ResultSet rs = stmt.executeQuery("select * from image_Shape_Feature");
		int a = 0;
		while (rs.next()) {
			a++;
			Original_fv[0][a]=rs.getDouble(2);
			Original_fv[1][a]=rs.getDouble(3);
			Original_fv[2][a]=rs.getDouble(8);
		}

		build_KD_tree(1,a,0,0);
		for(int i=1;i<=14;i++)
			DBUtil.node_id[i]=new int[level_size[i]+1];
		dep_to_node[1]=0;
		for(int i=1;i<=node_cnt;i++)
			if(Ls[i]!=-1) {
				dep_to_node[Ls[i]]=dep_to_node[Rs[i]]=dep_to_node[i]+1;
				DBUtil.node_id[dep_to_node[Ls[i]]][Ls_index[i]]=Ls[i];
				DBUtil.node_id[dep_to_node[Ls[i]]][Ls_index[i]+1]=Rs[i];
			}

		/*minimal_distance=1e12;
		minimal_id=-1;
		Tm=0;
		search_KD_tree(1);*/
	}
	public static Pair<Integer,Integer> build_KD_tree(int l, int r, int depth, int Fa_id) {
		if(l>r) return new Pair<>(-1,-1);
		Fa[++node_cnt]=Fa_id;
		level_size[depth]++;
		if(depth>mx_depth) mx_depth=depth;
		split_dim[node_cnt]=-1;
		if(l==r) {
			Ls[node_cnt]=Rs[node_cnt]=-1;
			Fv[0][node_cnt]=Original_fv[0][l];
			Fv[1][node_cnt]=Original_fv[1][l];
			Fv[2][node_cnt]=Original_fv[2][l];
			return new Pair<>(node_cnt, level_size[depth]);
		}
		int root_id=node_cnt;
		double mn_val=1e12;
		for(int dim=0;dim<3;dim++) {
			double average=0.0;
			for(int i=l;i<=r;i++)
				average+=Original_fv[dim][i];
			average/=(r-l+1);
			double val=0.0;
			for(int i=l;i<=r;i++)
				val+=(Original_fv[dim][i]-average)*(Original_fv[dim][i]-average);
			if(val<mn_val) {
				mn_val=val;
				split_dim[root_id]=dim;
			}
		}
		for(int i=l;i<=r;i++)
			tmp_fv[split_dim[root_id]][i]=Original_fv[split_dim[root_id]][i];
		Arrays.sort(tmp_fv[split_dim[root_id]],l,r+1);
		int L_point=l,R_point=r,mid=(l+r)>>1;
		double root_val=tmp_fv[split_dim[root_id]][mid];
		Fv[split_dim[root_id]][root_id]=root_val;
		for(int i=l;i<=r;i++) {
			int dir=(Original_fv[split_dim[root_id]][i]<=root_val)?1:2;
			if(dir==1) {
				tmp_fv[0][L_point]=Original_fv[0][i];
				tmp_fv[1][L_point]=Original_fv[1][i];
				tmp_fv[2][L_point]=Original_fv[2][i];
				L_point++;
			} else {
				tmp_fv[0][R_point]=Original_fv[0][i];
				tmp_fv[1][R_point]=Original_fv[1][i];
				tmp_fv[2][R_point]=Original_fv[2][i];
				R_point--;
			}
		}
		for(int i=l;i<=r;i++) {
			Original_fv[0][i]=tmp_fv[0][i];
			Original_fv[1][i]=tmp_fv[1][i];
			Original_fv[2][i]=tmp_fv[2][i];
		}
		Pair<Integer,Integer> L_info=build_KD_tree(l,mid,depth+1,root_id);
		Pair<Integer,Integer> R_info=build_KD_tree(mid+1,r,depth+1,root_id);
		Ls[root_id]=L_info.getKey();
		Rs[root_id]=R_info.getKey();
		Ls_index[root_id]=L_info.getValue();
		return new Pair<>(root_id,level_size[depth]);
	}

	public static int OT_Round;
	public static void final_test(int count)throws Exception {
		minimal_distance=1e12;
		minimal_id=-1;
		Fv_u=new double[3];
		ans_Fv=new double[3];
		//search_KD_tree_database(1);
		{
			long startTime=System.currentTimeMillis();
			OT_Round = 0;
			totalti=0;
			for(int i=0;i<1;i++) {
				DoOriginalImage();
				u_SPLIT_DIM=split_dim[1];
				u_FV=Fv[split_dim[1]][1];
				secure_search_KD_tree_database(1,0);
			}
			long endTime=System.currentTimeMillis();
			System.out.println("程序运行时间： "+(endTime-startTime)+"ms");
			System.out.println("totalti: "+totalti);
			System.out.println("OT_Round: "+OT_Round);
			System.out.println(ans_Fv[0]);
			System.out.println(ans_Fv[1]);
			System.out.println(ans_Fv[2]);
		}

		ResultSet=new float[count];
		ResultIDSet=new int[count];
		ResultPathSet=new String[count];
		DoComparison();
		mysort();
		ResultSet = ReturnSimilaritySet();
		for (int i = 0; i < OriginalImageShape.length; i++) {
			System.out.print(OriginalImageShape[i] + " | ");
		}
		System.out.println();
		for (int i = 0; i < 1; i++) {
			System.out.print("ResultSet[" + i + "]=" + ResultSet[i] + " , ");
			System.out.println("ResultIDSet[" + i + "]=" + ResultIDSet[i]);
			System.out.println("ResultPathSet[" + i + "]=" + ResultPathSet[i]);
		}
		System.out.println("相似度最大：" + ResultSet[0]);
	}

	public static void search_KD_tree(int u) {
		Tm++;

		if(split_dim[u]==-1) {
			double sum=0.0;
			for(int i=0;i<3;i++)
				sum+=(OriginalImageShape[i]-Fv[i][u])*(OriginalImageShape[i]-Fv[i][u]);
			double Dis=Math.sqrt(sum);
			if(Dis<minimal_distance) {
				minimal_id=u;
				minimal_distance=Dis;
			}
			return;
		}
		double Dis=OriginalImageShape[split_dim[u]]-Fv[split_dim[u]][u];
		if(Dis<0) {
			search_KD_tree(Ls[u]);
			if(Math.abs(Dis)<minimal_distance) search_KD_tree(Rs[u]);
		} else {
			search_KD_tree(Rs[u]);
			if(Math.abs(Dis)<minimal_distance) search_KD_tree(Ls[u]);
		}
	}

	public static int u_SPLIT_DIM;
	public static double u_FV;

	public static void Update_Inf_u(byte[] res) {
		u_SPLIT_DIM=ByteUtils.byteArrayToInt(res,0);
		u_FV=ByteUtils.bytes2Double(res,4);
	}

	private static int totalti=0;
	public static void Derandomization(int level,int desp) throws SQLException, NoSuchAlgorithmException, IOException {
		byte[] y = new byte[(DBUtil.node_id[level].length-1)*12];
		DBUtil.Set_Y(level,(desp+N-Selectionindexes[OT_Round])%N);
		DBUtil.Y_TransferSender();
		dr.readFully(y,0,(DBUtil.node_id[level].length-1)*12);
		totalti+=DBUtil.node_id[level].length-1;
		byte[] yy=new byte[12];
		System.arraycopy(y,desp*12,yy,0,12);
		byte[] hash=RandomOracles.G(12, t0jArray[OT_Round]);
		Update_Inf_u(xorByteArrays(hash,yy));
		OT_Round++;
	}

	public static void secure_search_KD_tree_database(int u, int level) throws SQLException, NoSuchAlgorithmException, IOException {
		int split_dim=u_SPLIT_DIM;
		if(split_dim==-1) {
			double sum=0.0;
			Fv_u[0]=Fv[0][u];
			Fv_u[1]=Fv[1][u];
			Fv_u[2]=Fv[2][u];
			for(int i=0;i<3;i++)
				sum+=(OriginalImageShape[i]-Fv_u[i])*(OriginalImageShape[i]-Fv_u[i]);
			double Dis=Math.sqrt(sum);
			if(Dis<minimal_distance) {
				minimal_id=u;
				minimal_distance=Dis;
				ans_Fv[0]=Fv_u[0];
				ans_Fv[1]=Fv_u[1];
				ans_Fv[2]=Fv_u[2];
			}
			return;
		}
		double Dis;
		Dis=OriginalImageShape[split_dim]-u_FV;
		int tmp_LS=Ls[u],tmp_RS=Rs[u],tmp_LS_INDEX=Ls_index[u];
		if(Dis<0) {
			Derandomization(level+1,tmp_LS_INDEX-1);
			secure_search_KD_tree_database(tmp_LS,level+1);
			if(Math.abs(Dis)<minimal_distance) {
				Derandomization(level+1,tmp_LS_INDEX);
				secure_search_KD_tree_database(tmp_RS,level+1);
			}
		} else {
			Derandomization(level+1,tmp_LS_INDEX);
			secure_search_KD_tree_database(tmp_RS,level+1);
			if(Math.abs(Dis)<minimal_distance) {
				Derandomization(level+1,tmp_LS_INDEX-1);
				secure_search_KD_tree_database(tmp_LS,level+1);
			}
		}
	}

	private static byte[] sArray;
	private static byte[] InputSender;//right now l=1.
	private static LinkedList<byte[]> InputSList=new LinkedList<byte[]>();

	public static void PreProcessing_Sender(int m,int n,int k,int l) {
		//Set InputList
		for(int j=0; j<m; j++) {
			InputSender=new byte[n*l/8];
			new Random().nextBytes(InputSender);
			InputSList.add(InputSender);
		}


		//Set sArray
		sArray=new byte[k/8];
		new Random().nextBytes(sArray);
	}


	void Preprocessing() {
		new Random().nextBytes(sArray);
	}

	private static BigInteger ddh_keys[]=null;
	private static BigInteger ddh_a=null;
	public static void BaseDDH_ithOT_GenerateKeys(int i) throws IOException {
		ddh_a=ElGamalEncryption.getRandom();
		ddh_keys=ElGamalEncryption.generatePublicAndFakeKeys(ddh_a);

		int bit=getBit(sArray, i);

		bw.write(ddh_keys[bit].toString());
		bw.newLine();
		bw.write(ddh_keys[1-bit].toString());
		bw.newLine();
		bw.flush();
	}

	private static byte[][] aArray=new byte[K][];
	public static void BaseDDH_ithOT_ReceiveMsg(int i,int k) throws IOException {
		BigInteger cipher0[]=new BigInteger[2];
		BigInteger cipher1[]=new BigInteger[2];
		cipher0[0]=new BigInteger(br.readLine());
		cipher0[1]=new BigInteger(br.readLine());
		cipher1[0]=new BigInteger(br.readLine());
		cipher1[1]=new BigInteger(br.readLine());

		int bit=getBit(sArray, i);
		BigInteger msg;
		if(bit==0) {
			msg=ElGamalEncryption.decrypt(cipher0[0], cipher0[1], ddh_a);
		} else {
			msg=ElGamalEncryption.decrypt(cipher1[0], cipher1[1], ddh_a);
		}
		aArray[i]= msg.toByteArray();
		if(aArray[i].length<k/8) {
			byte[] tmp = new byte[k/8-aArray[i].length];
			for(int j=0;j<k/8-aArray[i].length;j++)
				tmp[j]=0;
			byte[] tmp2=ByteUtils.byteMerger(tmp,aArray[i]);
			aArray[i] = tmp2;
		}
		else if (aArray[i][0] == 0&&aArray[i].length>k/8) {
			byte[] tmp = new byte[k/8];
			System.arraycopy(aArray[i], aArray[i].length-k/8, tmp, 0, tmp.length);
			aArray[i] = tmp;
		}
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
	private static LinkedList<byte[]> Y= new LinkedList<byte[]>();
	public static void setQ_QjArray_andY(int m, int n, int k, int l) throws NoSuchAlgorithmException {
		System.out.println();
		System.out.println("\t\tSet Q matrix and Y in IKNP's Sender begins now...");
		System.out.println();

		for (int i = 0; i < k; i++) {
			if (getBit(sArray, i) == 0)  {
				qArray[i] =xorByteArrays(RandomOracles.G(m/8,aArray[i]),wArray_Recieved[i]);
			}
			else if (getBit(sArray, i) == 1) {
				qArray[i] = RandomOracles.G(m/8,aArray[i]);
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

		for(int j=0;j<m;j++){
			byte[]and;
			byte[]xor;
			byte[] tmp=InputSList.get(j);

			for(int nn=0; nn<n;nn++) {
				and=andByteArrays(DBUtil.CodeWords.get(nn),sArray);
				xor=xorByteArrays(qjArray[j],and);
				byte[] H_val=RandomOracles.G(l/8,xor);
				byte[] y_nn=xorByteArrays(tmp,nn*l/8,l/8,H_val);
				Y.add(y_nn);
			}
		}

		System.out.println();
		System.out.println("\t\tSet Q matrix and Y in IKNP's Sender ends now...");
		System.out.println();
	}

	public static void Y_TransferSender(int n,int j,int nn) throws IOException {
		dw.write(Y.get(j*n+nn));
		dw.flush();
	}

	public static LinkedList<byte[]> CodeWords;
	static void Generate_OT() throws Exception {
		DBUtil.Generate_OT();
		CodeWords = CodeFunction.PRC(N);
		PreProcessing2_Receiver(M,N,K);
		DBUtil.PreProcessing2_Sender(M,N,K,L);
		OT_Extension(M,N,K,L);

	}

	private static int[] Selectionindexes=new int[M];
	private static byte[][] uArray=new byte[K][];
	private static byte[][] vArray=new byte[K][];
	private static byte[][] djArray=new byte[M][];
	private static byte[][] diArray=new  byte[K][];
	private static byte[][] t0Array=new  byte[K][];
	private static byte[][] t1Array=new  byte[K][];
	private static byte[][] t0jArray=new  byte[M][];
	public static void PreProcessing2_Receiver(int m, int n, int k) throws NoSuchAlgorithmException {
		//random choice r=(r1,r_2,...)
		for(int j=0; j<m; j++) {
			Random rn = new Random();
			Integer randomNum = rn.nextInt(n);
			Selectionindexes[j]=randomNum;
		}

		//Set matrix U and V
		for(int i=0; i < k; i++) {
			byte[] tmp=InputSList.get(i);
			uArray[i]=new byte[Kappa/8];
			vArray[i]=new byte[Kappa/8];
			System.arraycopy(tmp,0,uArray[i],0,Kappa/8);
			System.arraycopy(tmp,Kappa/8,vArray[i],0,Kappa/8);
		}

		//Set djArray: Code Matrix
		for (int j= 0; j<m; j++) {
			djArray[j]= CodeWords.get(Selectionindexes[j]);
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

	private static byte[][] wArray=new byte[K][];
	public static void OT_Extension(int m,int n,int k,int l) throws IOException, NoSuchAlgorithmException {
		//Set wArray:t0 xor G(u)
		for (int i = 0; i < k; i++) {
			wArray[i]=(xorByteArrays(RandomOracles.G(m/8,uArray[i]),t0Array[i]));
		}

		{
			System.out.println();
			System.out.println("\t\tSending wArray in OPRF's Receiver begins now...");
			System.out.println();

			for (int i = 0; i < k; i++) {
				dw.writeInt(wArray[i].length);
				dw.write(wArray[i]);
				dw.flush();
				DBUtil.W_TransferReceiver(i);
			}

			System.out.println();
			System.out.println("\t\tSending wArray OPRF's Receiver ends now...");
			System.out.println();
		}

		DBUtil.setQ_QjArray(m,n,k,l);
	}

	/*
	 * 返回路径集
	 */
	public String[] ReturnPathSet(String OriginalPath) throws Exception {
		int count=new GetRecordCount().getCount();
		ResultSet=new float[count];
		ResultIDSet=new int[count];
		ResultPathSet=new String[count];
		DoOriginalImage(OriginalPath);
		DoComparison();
		mysort();
		return ResultPathSet;
	}
	/*
	 * 返回相似度集
	 */
	public static float[] ReturnSimilaritySet() {
		for (int i = 0; i < ResultSet.length; i++) {
			ResultSet[i] = ResultSet[i] / ResultSet[ResultSet.length - 1];
		}
		for (int i = 0; i < ResultSet.length; i++) {
			ResultSet[i] = (float) (1.0 - ResultSet[i]);
		}
		return ResultSet;
	}

	

	/*
	 * 实现地址、ID、距离差同步排序
	 */
	public static void mysort() {
		float temp;
		int tempid;
		String tempPath;
		for (int i = 0; i < ResultSet.length; i++) {
			for (int j = i; j < ResultSet.length; j++) {
				if (ResultSet[i] > ResultSet[j]) {
					temp = ResultSet[i];
					ResultSet[i] = ResultSet[j];
					ResultSet[j] = temp;
					tempid = ResultIDSet[i];
					ResultIDSet[i] = ResultIDSet[j];
					ResultIDSet[j] = tempid;
					tempPath = ResultPathSet[i];
					ResultPathSet[i] = ResultPathSet[j];
					ResultPathSet[j] = tempPath;
				}
			}
		}
	}

	/*
	 * 获取地址结果集
	 */
	public static void SetResultPathSet() throws SQLException {
		Connection con = DBUtil.getConnection();
		Statement stmt = con.createStatement();
		ResultSet rs = stmt.executeQuery("select imagePath from image_Library");
		int i = 0;
		while (rs.next()) {
			ResultPathSet[i++] = rs.getString("imagePath");
		}
	}

	/*
	 * 求距离逼存到ResultSet数组里
	 */
	public static float SeekingDistance() {
		float ret = 0;
		double sum = 0;
		for (int i = 0; i < OriginalImageShape.length; i++) {
			sum = sum + Math.pow((OriginalImageShape[i] - ComparisonImageShape[i]), 2);
		}
		ret = (float) Math.sqrt(sum / OriginalImageShape.length);
		return ret;
	}

	/*
	 * 原图特征值和数据库中的对比
	 */
	public static void DoComparison() throws Exception {
		Connection con = DBUtil.getConnection();
		Statement stmt = con.createStatement();
		ResultSet rs = stmt.executeQuery("select * from image_Shape_Feature");
		int a = 0;

		while (rs.next()) {
			for (int i = 0; i < ComparisonImageShape.length - 1; i++) {
				ComparisonImageShape[i] = rs.getDouble(i + 2);
			}
			ComparisonImageShape[ComparisonImageShape.length - 1] = rs.getDouble(8);
			ResultSet[a] = SeekingDistance();
			ResultIDSet[a] = ++a;
		}
		SetResultPathSet();
	}

	/*
	 * 求原图特征值并存到OriginalImageShape数组里
	 */
	public static void DoOriginalImage() throws IOException, SQLException {
		System.out.print("请输入原始图的地址：");
		BufferedReader br2 = new BufferedReader(new InputStreamReader(System.in));
		String OriginalImagePath = br2.readLine();
		//System.out.println(OriginalImagePath);
		OriginalImagePath = OriginalImagePath.replace("\\", "\\/");
		double shape[] = new double[7];
		shape = GetImageShapeFeature.getShapeFeature(OriginalImagePath);
		OriginalImageShape[0] = shape[0];
		OriginalImageShape[1] = shape[1];
		OriginalImageShape[2] = shape[6];
	}

	public static void DoOriginalImage(String OriginalImagePath) throws IOException {
		OriginalImagePath = OriginalImagePath.replace("\\", "\\/");
		double shape[] = new double[7];
		shape = GetImageShapeFeature.getShapeFeature(OriginalImagePath);
		OriginalImageShape[0] = shape[0];
		OriginalImageShape[1] = shape[1];
		OriginalImageShape[2] = shape[6];
	}



	public static byte[] xorByteArrays(byte[] a, int from, int len, byte[] b) {

		if (!(len == b.length)) {
			System.out.println("Lengths NOT equal");
			System.exit(3);
		}
		byte[] result = new byte[len];

		for (int i = 0; i < result.length; i++) {
			result[i] = (byte) (((int) a[i+from]) ^ ((int) b[i]));
		}
		return result;
	}


	// Method to get the and result between two byte arrays

	public static byte[] andByteArrays(byte[] a, byte[] b) {

		byte[] result = new byte[Math.min(a.length, b.length)];

		if (!(a.length == b.length)) {
			System.out.println("Lengths NOT equal");
			System.exit(3);
		}

		for (int i = 0; i < result.length; i++) {
			result[i] = (byte) (((int) a[i]) & ((int) b[i]));
		}
		return result;
	}


	// Method to convert an int array into a byte array

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


	//  Method to convert the byte array into binary strings

	public String toBinary( byte[] bytes )
	{
		StringBuilder sb = new StringBuilder(bytes.length * Byte.SIZE);
		for( int i = 0; i < Byte.SIZE * bytes.length; i++ )
			sb.append((bytes[i / Byte.SIZE] << i % Byte.SIZE & 0x80) == 0 ? '0' : '1');
		return sb.toString();
	}

	// Method to read the value of a bit (0 or 1) of a byte array

	private static int getBit(byte[] data, int pos) {
		int posByte = pos/8;
		int posBit = pos%8;
		byte valByte = data[posByte];
		int valInt = valByte>>(8-(posBit+1)) & 0x0001;
		return valInt;
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
		/*

		System.out.print(OriginalImageShape[0]);
		System.out.print(' ');
		System.out.print(OriginalImageShape[1]);
		System.out.print(' ');
		System.out.println(OriginalImageShape[2]);

		System.out.print(Fv[0][minimal_id]);
		System.out.print(' ');
		System.out.print(Fv[1][minimal_id]);
		System.out.print(' ');
		System.out.println(Fv[2][minimal_id]);*/

		/*PrintStream output_sql=new PrintStream("/Users/yangheng/Desktop/CBIR-master/shape_kd_tree.txt");
		for(int i=1;i<=node_cnt;i++) {
			output_sql.print("INSERT INTO `shape_kd_tree` VALUES ('");
			output_sql.print(i);
			output_sql.print("', '");
			output_sql.print(Ls[i]);
			output_sql.print("', '");
			output_sql.print(Rs[i]);
			output_sql.print("', '");
			output_sql.print(Fa[i]);
			output_sql.print("', '");
			output_sql.print(split_dim[i]);
			output_sql.print("', '");
			output_sql.print(Ls_index[i]);
			output_sql.print("', '");
			output_sql.print(Fv[0][i]);
			output_sql.print("', '");
			output_sql.print(Fv[1][i]);
			output_sql.print("', '");
			output_sql.print(Fv[2][i]);
			output_sql.println("');");
		}*/

/*public static void search_KD_tree_database(int u) throws SQLException {
		Connection con = DBUtil.getConnection();
		Statement stmt = con.createStatement();
		ResultSet rs = stmt.executeQuery("select * from shape_kd_tree");
		rs.absolute(u);
		int split_dim=rs.getInt("SPLIT_DIM");
		if(split_dim==-1) {
			double sum=0.0;
			Fv_u[0]=rs.getDouble("FV0");
			Fv_u[1]=rs.getDouble("FV1");
			Fv_u[2]=rs.getDouble("FV2");
			for(int i=0;i<3;i++)
				sum+=(OriginalImageShape[i]-Fv_u[i])*(OriginalImageShape[i]-Fv_u[i]);
			double Dis=Math.sqrt(sum);
			if(Dis<minimal_distance) {
				minimal_id=u;
				minimal_distance=Dis;
				ans_Fv[0]=Fv_u[0];
				ans_Fv[1]=Fv_u[1];
				ans_Fv[2]=Fv_u[2];
			}
			return;
		}
		double Dis=OriginalImageShape[split_dim]-rs.getDouble("Fv"+(char)('0'+split_dim));
		if(Dis<0) {
			int Ls=rs.getInt("LS");
			search_KD_tree_database(Ls);
			if(Math.abs(Dis)<minimal_distance) {
				int Rs=rs.getInt("RS");
				search_KD_tree_database(Rs);
			}
		} else {
			int Rs=rs.getInt("RS");
			search_KD_tree_database(Rs);
			if(Math.abs(Dis)<minimal_distance) {
				int Ls=rs.getInt("LS");
				search_KD_tree_database(Ls);
			}
		}
	}*/

