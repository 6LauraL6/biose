'use client'
// Navbar.tsx
// Navbar.tsx
import { useState } from "react";
import Link from "next/link";
import Logo from "./Logo";

const Navbar = () => {
  const [isMenuOpen, setMenuOpen] = useState(false);

  const toggleMenu = () => {
    setMenuOpen(!isMenuOpen);
  };

  return (
    <>
      <div className="w-full h-20 bg-indigo-300 sticky top-0">
        <div className="container mx-auto px-4 h-full">
          <div className="flex justify-between items-center h-full">
            <div className="md:hidden">
              {/* Utilitzem l'atribut onClick directament */}
              <button className="text-white" onClick={toggleMenu}>
                <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24"
                  stroke-linecap="round" stroke-linejoin="round" stroke-width="2">
                  <path d="M4 6h16M4 12h16m-7 6h7"></path>
                </svg>
              </button>
              <ul className={`md:hidden ${isMenuOpen ? '' : 'hidden'} gap-x-6 text-white absolute top-20 left-0 right-0 bg-black p-2 rounded shadow-md`}>
                <li className="mb-2">
                  <Link href="/">
                    <a className="text-white hover:text-red-500 block">Main</a>
                  </Link>
                </li>
                <li className="mb-2">
                  <Link href="/act1">
                    <a className="text-white hover:text-red-500 block">Act1</a>
                  </Link>
                </li>
                <li>
                  <Link href="/act2">
                    <a className="text-white hover:text-red-500 block">Act2</a>
                  </Link>
                </li>
              </ul>
            </div>
            <ul className="hidden md:flex gap-x-6 text-white">
              <Logo />
              <li>
                <Link href="/">
                  <a>Main</a>
                </Link>
              </li>
              <li>
                <Link href="/act1">
                  <a>Act1</a>
                </Link>
              </li>
              <li>
                <Link href="/act2">
                  <a>Act2</a>
                </Link>
              </li>
            </ul>
          </div>
        </div>
      </div>
    </>
  );
};

export default Navbar;
