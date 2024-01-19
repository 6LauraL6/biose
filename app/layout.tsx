import type { Metadata } from 'next'
import { Inter } from 'next/font/google'
import Navigation from './components/navigation/navbar'
import './globals.css'
import './locals.css'

const inter = Inter({ subsets: ['latin'] })

export const metadata: Metadata = {
  title: 'Bio - Sequence',
  description: '',
}

export default function RootLayout({
  children,
}: {
  children: React.ReactNode
}) {
  return (
    <html>
      <head>
        <title>Bio Sequence projects.</title>
        <meta name="viewport" content="width=device-width, initial-scale=1" />
        <link rel="stylesheet" href=" https://cdn.jsdelivr.net/npm/bootstrap@5.3.1/dist/css/bootstrap.min.css" />
      </head>
      <body className={inter.className}>
        <Navigation />
        {children}
        <footer className="container">
          <p className="text-center">
            <a href="https://creativecommons.org/licenses/by-nc/4.0/deed.ca">CC BY-NC 4.0 Deed</a>  —  David de Mingo, Miquel Amorós
          </p>
        </footer>
      </body>
    </html>
  )
}
