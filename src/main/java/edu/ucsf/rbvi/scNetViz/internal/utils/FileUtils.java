package edu.ucsf.rbvi.scNetViz.internal.utils;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;

import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;

public class FileUtils {
	public static void skipHeader(List<String[]> table) {
		if ((table.size() > 0) && (table.get(0).length > 1) && 
		    !ModelUtils.isInteger(table.get(0)[0])) {
			// It is -- skip it
			System.out.println("Skipping header");
			table.remove(0);
		}
		return;
	}

	public static boolean isTar(String fileName) {
		if (fileName.endsWith(".tar") || fileName.endsWith(".tgz") || fileName.endsWith(".tar.gz"))
			return true;
		return false;
	}

	public static boolean isZip(String fileName) {
		if (fileName.endsWith(".zip"))
			return true;
		return false;
	}

	public static boolean isGzip(String fileName) {
		if (fileName.endsWith(".gz"))
			return true;
		return false;
	}

	public static InputStream getGzipStream(InputStream stream) throws IOException {
		// return new GZIPInputStream(stream);
		return new GzipCompressorInputStream(stream);
	}

	public static ZipInputStream getZipInputStream(File mtxFile) throws IOException {
		FileInputStream inputStream = new FileInputStream(mtxFile);
		ZipInputStream zipStream = new ZipInputStream(inputStream);
		return zipStream;
	}

	public static TarArchiveInputStream getTarInputStream(File mtxFile) throws IOException {
		InputStream inputStream = new FileInputStream(mtxFile);
		if (mtxFile.getName().endsWith(".gz") || mtxFile.getName().endsWith(".tgz")) {
			System.out.println("Getting a gzip input stream");
			inputStream = new GzipCompressorInputStream(inputStream);
		}
		TarArchiveInputStream tarStream = new TarArchiveInputStream(inputStream);
		return tarStream;
	}

	public static String baseName(String name) {
		if (name.contains(File.separator)) {
			name = name.substring(name.lastIndexOf(File.separator)+1);
		}
		if (name.endsWith(".tar.gz"))
			return name.substring(0,name.length()-7);
		if (name.endsWith(".gz"))
			return name.substring(0,name.length()-3);
		if (name.endsWith(".zip") || name.endsWith(".tgz"))
			return name.substring(0,name.length()-4);
		return name;
	}

}
