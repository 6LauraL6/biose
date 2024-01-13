'use client';
import seq from "bionode-seq";
import { useEffect, useState } from "react";

export default function Board() {

  interface GenbankInfo {
    title: string;
    accession: string;
    ncbi_link: string;
    latest_reference: string;
    num_features: number;
    cds_info: { type: string; location: string }[];
    sequence: string;
    gc: number;
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

      // Registra la resposta per a fins de depuració
      console.log('Response from server:', data);
  
      // Intenta aplanar l'array, assumeixent que és un array de genbank_info
      const allGenbankInfo = Array.isArray(data) ? data.flat() : [];
      setGenbankInfo(allGenbankInfo);
    } catch (error) {
      console.error("Error fetching Genbank info:", error);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    fetchGenbankInfo();
  }, []);

  return (
    <main className="p-6">
  <h1 className="text-indigo-700 text-xl font-bold">BioActivitat 4 - Fitxers de seqüències.</h1>
  <p>⚠ Per ara només tenim allotjats fitxers Unigenbank de com a molt 1 MB d'espai. ⚠</p>
  <div className="bg-green shadow-md rounded px-8 pt-6 pb-8 mb-4 grid grid-cols-1 lg:grid-cols-2">
    {loading && (
      <div className="flex items-center justify-center mt-4">
        <div className="border-t-8 border-gray-300 border-solid border-b-8 border-gray-300 rounded-full h-16 w-16 animate-spin"></div>
        <p className="ml-4">Carregant dades...</p>
      </div>
    )}
    {error && <p>{error}</p>}
    {genbankInfo.map((info, index) => (
      <div key={index} className="p-2">
        <h2 className="text-l font-bold">{info.title}</h2>
        <p>Accession ID: {info.accession}</p>
        <p className="text-indigo-700 no-underline hover:underline">🔗 Enllaç NCBI: 
          <a href={info.ncbi_link} target="_blank" rel="noopener noreferrer">{info.accession}</a>
        </p>
        <p>🖧 Referència Més Recent: {info.latest_reference}</p>
        <p>🧫 Número de Features: {info.num_features}</p>
        <p>🧬 Seqüència (30 primeres): {info.sequence.slice(0, 30)}</p>
        <p>🧬 %GC: {info.gc}</p>
        {/* Informació CDS */}
        {info.cds_info && info.cds_info.slice(0, 5).length > 0 && (
          <p>🔬 Informació primers CDS:
            <ul>
              {info.cds_info.slice(0, 5).map((cds, cdsIndex) => (
                <li key={cdsIndex}>{`Type: ${cds.type}, Location: ${cds.location}`}</li>
              ))}
            </ul>
          </p>
        )}
        <hr />
      </div>
    ))}
  </div>
</main>

  );
}