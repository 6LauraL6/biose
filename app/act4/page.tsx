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
  <h1 className="text-indigo-700 text-xl font-bold">BioActivitat 4 - Fitxer Genbank.</h1>
  <p>⚠ Per ara només tenim allotjats fitxers Unigenbank de com a molt 5 MB d&apos; espai. ⚠</p>
  <div className="bg-green shadow-md rounded px-8 pt-6 pb-8 mb-4 grid grid-cols-1 lg:grid-cols-2">
    {loading && (
      <div className="flex items-center justify-center mt-4">
        <div className="border-t-8 border-gray-300 border-solid border-b-8 border-gray-300 rounded-full h-16 w-16 animate-spin"></div>
        <p className="ml-4">Carregant dades...</p>
      </div>
    )}
    {error && <p>{error}</p>}
    
      <div className="p-2">
        <h2 className="text-l font-bold">Dades del Genbank</h2>
        
      </div>
    
  </div>
</main>
  );
}