'use client';
import seq from "bionode-seq";
import { useEffect, useState } from "react";

export default function Board() {

  // Necessitem aquest objecte.
  interface GenbankInfo {
    title: string;
    accession: string;
    ncbi_link: string;
    latest_reference: string;
    num_features: number;
    seq_first30: string;
    cds_info: { type: string; location: string }[];
    repeated_sequence: string;
  }

  const [loading, setLoading] = useState(false);
  const [genbankInfo, setGenbankInfo] = useState<GenbankInfo[]>([]);
  const [error, setError] = useState(null);

  // Funció per obtenir la informació del genbank
  const fetchGenbankInfo = async () => {
    setLoading(true);
  
    try {
      const response = await fetch('/api/genbank-list');
      const data = await response.json();
      const allGenbankInfo = data.flat();
      setGenbankInfo(allGenbankInfo);
    } catch (error) {
      console.error("Error fetching Genbank info:", error);
      // setError("Error fetching Genbank info");
    } finally {
      setLoading(false);
    }
  };


  useEffect(() => {
    fetchGenbankInfo();
  }, [])

  return (
    <main className="p-6">
      <h1 className="text-indigo-700 text-xl font-bold">BioActivitat 4 - Fitxers de seqüències.</h1>
      <div className="bg-green shadow-md rounded px-8 pt-6 pb-8 mb-4 basis-1/2">
        {/* <button className="bg-blue-500 text-white py-2 px-4 rounded-full hover:bg-blue-600 focus:outline-none focus:shadow-outline-blue"
          onClick={fetchGenbankInfo}>Obtenir Informació Genbank
        </button> */}
        {loading && (
          <div className="flex items-center justify-center mt-4">
            <div className="border-t-8 border-gray-300 border-solid border-b-8 border-gray-300 rounded-full h-16 w-16 animate-spin"></div>
            <p className="ml-4">Carregant dades...</p>
          </div>
        )}
        {error && <p>{error}</p>}
        {genbankInfo && genbankInfo.map((info, index) => (
          <div key={index}>
            <h2 className="text-l font-bold">{info.title}</h2>
            <p className="text-l font-bold">Accession ID: {info.accession}</p>
            <p>Enllaç NCBI: <a className="no-underline hover:underline text-indigo-700" href={info.ncbi_link} target="_blank" rel="noopener noreferrer">Link NCBI: {info.accession}</a></p>
            <p>Referència Més Recent: {info.latest_reference}</p>
            <p>Número de Funcionalitats: {info.num_features}</p>
            {/* Informació CDS (amb comprovació si existeix) */}
            {info.cds_info && (
              <div>
                  <p>Informació CDS:</p>
                  <ul>
                    {info.cds_info.map((cds, cdsIndex) => (
                        <li key={cdsIndex}>{`Type: ${cds.type}, Location: ${cds.location}`}</li>
                    ))}
                  </ul>
              </div>
            )}
            {/* Informació Seqüència Repetida */}
            {info.repeated_sequence && (
              <p>
                Seqüència Repetida: {info.repeated_sequence}
              </p>
            )}
            <hr/>
          </div>
      ))}
      </div>
    </main>
  );
}